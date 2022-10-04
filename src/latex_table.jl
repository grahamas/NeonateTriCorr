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

function prettify_df(df, row_fn)
    DataFrame(row_fn.(eachrow(df)))
end

function fixed_reduction_rev_set_excl_row(r)
    (
        target=r.target,
        standardization = r.standardization,
        signal = r.signal,
        AUC = r.auc,
        TPR = r.tpr57,
        FPHr = r.fphr80
    )
end

using LaTeXStrings

function latexify_sensitivity(df::DataFrame)
    """
Target & Standardization & Signal & AUC & TPR (5.7 FP/Hr) & FP/Hr (80\\% TPR)\\\\
    \\hline
    """ *
    join(latexify_sensitivity.(eachrow(df)), "\\\\\n")
end

round_if_number(x) = x
round_if_number(x::Number) = round(x, sigdigits=2) 

function latexify_sensitivity(row::DataFrameRow)
    row_contents = join(round_if_number.(getproperty.(Ref(row), ["target", "standardization", "signal", "auc", "tpr57", "fphr80"])), "&")
    return row_contents
end
