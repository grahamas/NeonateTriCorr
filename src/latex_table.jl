function name_info_flow_row(row)
    (spike_motif=row[1],  n1_num=row[2], t1_num=row[3], n2_num=row[4],  t2_num=row[5], motif_class_num=row[6])
end

sign_str(x) = if x > 0
    "+"
elseif x < 0
    "-"
else
    "0"
end
abs_rel_sym(x1, x2) = if abs(x1) < abs(x2)
    "<"
elseif abs(x1) > abs(x2)
    ">"
else
    "="
end

function conditions_from_lag_nums(n1, t1, n2, t2)
    conditions = filter(x -> x != "", [
        if sign_str(n1) == sign_str(n2) && n1 != 0
            "\$\\abs{n_1} " * abs_rel_sym(n1, n2) * " \\abs{n_2}\$"
        else
            ""
        end, if sign_str(t1) == sign_str(t2) && t1 != 0
            "\$\\abs{t_1} " * abs_rel_sym(t1, t2) * " \\abs{t_2}\$"
        else
            ""
        end])
    conditions_str = if length(conditions) > 1
        "\\adjustbox{raise=-0.25\\height}{\\shortstack{" * join(conditions, "\\\\") * "}}"
    elseif length(conditions) == 0
        ""
    else
        only(conditions)
    end
    return conditions_str
end

function prettify_df(df)
    DataFrame(prettify_row.(eachrow(df)))
end

function prettify_row(row)
    (
        spike_motif = row.spike_motif,
        n1 = sign_str(row.n1_num),
        t1 = sign_str(row.t1_num),
        n2 = sign_str(row.n2_num),
        t2 = sign_str(row.t2_num),
        motif_class = roman_encode(row.motif_class_num),
        conditions = conditions_from_lag_nums(row.n1_num, row.t1_num, row.n2_num, row.t2_num)
    )
end

using LaTeXStrings

function latexify(df::DataFrame)
    """
Lag Motif & \$n_1\$ & \$t_1\$ & \$n_2\$ & \$t_2\$ & Constraints & Information Flow Pattern & Configuration \\\\
    \\hline
    """ *
    join(latexify.(eachrow(df)), "\\\\\n")
end

function graphics_scaling_string(n_same_neurons, n_same_times)
    if n_same_neurons == 3 && n_same_times == 3
        return "height=1.5cm"
    else
        return "scale=1.0"
    end
end

function latexify(row::DataFrameRow)
    row_contents = join(getproperty.(Ref(row), ["spike_motif", "n1", "t1", "n2", "t2", "conditions", "motif_class"]), "&")
    row_contents *= "& \\adjustbox{padding=0em 0.5em 0em 0.5em, raise=-0.33\\height}{\\includegraphics{motif_figs/$(row.spike_motif).png}}"
    return row_contents
end

function remap_old_filenames(srcdir, dstdir)
    old_classes = generate_all_lag_motif_classes_OLD()
    new_classes = generate_all_lag_motif_classes()
    @show old_classes
    @show new_classes
    ntnt2newname_mapping = Dict(
        tup[2:end] => tup[1] for tup in new_classes
    )
    @show ntnt2newname_mapping
    oldname2newname_mapping = Dict(
        tup[1] => ntnt2newname_mapping[(tup[2],tup[4],tup[3],tup[5])] for tup in old_classes
    )
    @show oldname2newname_mapping
    for (src, dst) in pairs(oldname2newname_mapping)
        cp(joinpath(srcdir, "$(src).png"), joinpath(dstdir, "$(dst).png"))
    end
end
