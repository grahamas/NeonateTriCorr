quote
    begin
        #= DIR/loopvec_test.jl:24 =#
        nothing
    end
    begin
        #= DIR/loopvec_test.jl:22 =#
        nothing
    end
    begin
        #= DIR/loopvec_test.jl:21 =#
        nothing
    end
    begin
        #= DIR/loopvec_test.jl:20 =#
        nothing
    end
    begin
        #= DIR/loopvec_test.jl:19 =#
        nothing
    end
    begin
        #= DIR/loopvec_test.jl:19 =#
        nothing
    end
    begin
        #= DIR/loopvec_test.jl:18 =#
        nothing
    end
    begin
        #= DIR/loopvec_test.jl:17 =#
        nothing
    end
    begin
        #= DIR/loopvec_test.jl:13 =#
        nothing
    end
    var"###looprangea1i###1###" = LoopVectorization.canonicalize_range(begin
                $(Expr(:inbounds, true))
                local var"#2187#val" = eachindex(range_a)
                $(Expr(:inbounds, :pop))
                var"#2187#val"
            end)
    var"###looplena1i###2###" = ArrayInterface.static_length(var"###looprangea1i###1###")
    var"###a1i_loop_lower_bound###3###" = LoopVectorization.maybestaticfirst(var"###looprangea1i###1###")
    var"###a1i_loop_upper_bound###4###" = LoopVectorization.maybestaticlast(var"###looprangea1i###1###")
    var"###looprangea2i###5###" = LoopVectorization.canonicalize_range(begin
                $(Expr(:inbounds, true))
                local var"#2188#val" = eachindex(range_a)
                $(Expr(:inbounds, :pop))
                var"#2188#val"
            end)
    var"###looplena2i###6###" = ArrayInterface.static_length(var"###looprangea2i###5###")
    var"###a2i_loop_lower_bound###7###" = LoopVectorization.maybestaticfirst(var"###looprangea2i###5###")
    var"###a2i_loop_upper_bound###8###" = LoopVectorization.maybestaticlast(var"###looprangea2i###5###")
    var"###looprangeb1i###9###" = LoopVectorization.canonicalize_range(begin
                $(Expr(:inbounds, true))
                local var"#2189#val" = eachindex(range_b)
                $(Expr(:inbounds, :pop))
                var"#2189#val"
            end)
    var"###looplenb1i###10###" = ArrayInterface.static_length(var"###looprangeb1i###9###")
    var"###b1i_loop_lower_bound###11###" = LoopVectorization.maybestaticfirst(var"###looprangeb1i###9###")
    var"###b1i_loop_upper_bound###12###" = LoopVectorization.maybestaticlast(var"###looprangeb1i###9###")
    var"###looprangeb2i###13###" = LoopVectorization.canonicalize_range(begin
                $(Expr(:inbounds, true))
                local var"#2190#val" = eachindex(range_b)
                $(Expr(:inbounds, :pop))
                var"#2190#val"
            end)
    var"###looplenb2i###14###" = ArrayInterface.static_length(var"###looprangeb2i###13###")
    var"###b2i_loop_lower_bound###15###" = LoopVectorization.maybestaticfirst(var"###looprangeb2i###13###")
    var"###b2i_loop_upper_bound###16###" = LoopVectorization.maybestaticlast(var"###looprangeb2i###13###")
    var"###looprangei_a###18###" = LoopVectorization.canonicalize_range(begin
                $(Expr(:inbounds, true))
                local var"#2191#val" = padded_axis_a
                $(Expr(:inbounds, :pop))
                var"#2191#val"
            end)
    var"###loopleni_a###19###" = ArrayInterface.static_length(var"###looprangei_a###18###")
    var"###i_a_loop_lower_bound###20###" = LoopVectorization.maybestaticfirst(var"###looprangei_a###18###")
    var"###i_a_loop_upper_bound###21###" = LoopVectorization.maybestaticlast(var"###looprangei_a###18###")
    var"###i_a_loop_step###22###" = LoopVectorization.static_step(var"###looprangei_a###18###")
    var"###looprangei_b###23###" = LoopVectorization.canonicalize_range(begin
                $(Expr(:inbounds, true))
                local var"#2192#val" = padded_axis_b
                $(Expr(:inbounds, :pop))
                var"#2192#val"
            end)
    var"###loopleni_b###24###" = ArrayInterface.static_length(var"###looprangei_b###23###")
    var"###i_b_loop_lower_bound###25###" = LoopVectorization.maybestaticfirst(var"###looprangei_b###23###")
    var"###i_b_loop_upper_bound###26###" = LoopVectorization.maybestaticlast(var"###looprangei_b###23###")
    var"###i_b_loop_step###27###" = LoopVectorization.static_step(var"###looprangei_b###23###")
    if LoopVectorization.check_args(range_a, range_b, src, target_arr)
        (var"##vptr##_range_a", var"#range_a#preserve#buffer#") = LoopVectorization.stridedpointer_preserve(range_a)
        (var"##vptr##_range_b", var"#range_b#preserve#buffer#") = LoopVectorization.stridedpointer_preserve(range_b)
        (var"##vptr##_src", var"#src#preserve#buffer#") = LoopVectorization.stridedpointer_preserve(src)
        (var"##vptr##_target_arr", var"#target_arr#preserve#buffer#") = LoopVectorization.stridedpointer_preserve(target_arr)
        var"####grouped#strided#pointer####39###" = Core.getfield(LoopVectorization.grouped_strided_pointer((LoopVectorization.densewrapper(LoopVectorization.gespf1(var"##vptr##_range_a", (LoopVectorization.similardims(var"###looprangea1i###1###", VectorizationBase.NullStep()),)), range_a), LoopVectorization.densewrapper(LoopVectorization.gespf1(var"##vptr##_range_b", (LoopVectorization.similardims(var"###looprangeb1i###9###", VectorizationBase.NullStep()),)), range_b), LoopVectorization.densewrapper(LoopVectorization.gespf1(var"##vptr##_src", (LoopVectorization.similardims(var"###looprangei_a###18###", VectorizationBase.NullStep()), LoopVectorization.similardims(var"###looprangei_b###23###", VectorizationBase.NullStep()))), src), LoopVectorization.densewrapper(LoopVectorization.gespf1(var"##vptr##_target_arr", (LoopVectorization.similardims(var"###looprangeb1i###9###", VectorizationBase.NullStep()), LoopVectorization.similardims(var"###looprangeb2i###13###", VectorizationBase.NullStep()))), target_arr)), Val{()}()), 1)
        $(Expr(:gc_preserve, quote
    var"##vargsym#1186" = ((var"###looprangea1i###1###", var"###looprangea2i###5###", var"###looprangeb1i###9###", var"###looprangeb2i###13###", var"###looprangei_a###18###", var"###looprangei_b###23###"), (var"####grouped#strided#pointer####39###",))
    var"##Tloopeltype##" = LoopVectorization.promote_type(LoopVectorization.eltype(range_a), LoopVectorization.eltype(range_b), LoopVectorization.eltype(src), LoopVectorization.eltype(target_arr))
    var"##Wvecwidth##" = LoopVectorization.pick_vector_width(var"##Tloopeltype##")
    LoopVectorization._turbo_!(LoopVectorization.avx_config_val(Val{(false, 0, 0, 0, false, 0x0000000000000001)}(), var"##Wvecwidth##"), Val{(:LoopVectorization, :getindex, LoopVectorization.OperationStruct(0x00000000000000000000000000000001, 0x00000000000000000000000000000000, 0x00000000000000000000000000000000, 0x00000000000000000000000000000000, 0x00000000000000000000000000000000, LoopVectorization.memload, 0x0001, 0x01), :LoopVectorization, :getindex, LoopVectorization.OperationStruct(0x00000000000000000000000000000002, 0x00000000000000000000000000000000, 0x00000000000000000000000000000000, 0x00000000000000000000000000000000, 0x00000000000000000000000000000000, LoopVectorization.memload, 0x0002, 0x02), :LoopVectorization, :getindex, LoopVectorization.OperationStruct(0x00000000000000000000000000000003, 0x00000000000000000000000000000000, 0x00000000000000000000000000000000, 0x00000000000000000000000000000000, 0x00000000000000000000000000000000, LoopVectorization.memload, 0x0003, 0x03), :LoopVectorization, :getindex, LoopVectorization.OperationStruct(0x00000000000000000000000000000004, 0x00000000000000000000000000000000, 0x00000000000000000000000000000000, 0x00000000000000000000000000000000, 0x00000000000000000000000000000000, LoopVectorization.memload, 0x0004, 0x04), :numericconstant, Symbol("###0###17###"), LoopVectorization.OperationStruct(0x00000000000000000000000000001234, 0x00000000000000000000000000000000, 0x00000000000000000000000000000056, 0x00000000000000000000000000000000, 0x00000000000000000000000000000000, LoopVectorization.constant, 0x0005, 0x00), :LoopVectorization, :getindex, LoopVectorization.OperationStruct(0x00000000000000000000000000000056, 0x00000000000000000000000000000000, 0x00000000000000000000000000000000, 0x00000000000000000000000000000000, 0x00000000000000000000000000000000, LoopVectorization.memload, 0x0006, 0x05), :i_a, :i_a, LoopVectorization.OperationStruct(0x00000000000000000000000000000005, 0x00000000000000000000000000000000, 0x00000000000000000000000000000000, 0x00000000000000000000000000000000, 0x00000000000000000000000000000000, LoopVectorization.loopvalue, 0x0007, 0x00), :LoopVectorization, :+, LoopVectorization.OperationStruct(0x00000000000000000000000000000015, 0x00000000000000000000000000000000, 0x00000000000000000000000000000000, 0x00000000000000000000000000010007, 0x00000000000000000000000000000000, LoopVectorization.compute, 0x0008, 0x00), :i_b, :i_b, LoopVectorization.OperationStruct(0x00000000000000000000000000000006, 0x00000000000000000000000000000000, 0x00000000000000000000000000000000, 0x00000000000000000000000000000000, 0x00000000000000000000000000000000, LoopVectorization.loopvalue, 0x0009, 0x00), :LoopVectorization, :+, LoopVectorization.OperationStruct(0x00000000000000000000000000000036, 0x00000000000000000000000000000000, 0x00000000000000000000000000000000, 0x00000000000000000000000000030009, 0x00000000000000000000000000000000, LoopVectorization.compute, 0x000a, 0x00), :LoopVectorization, :getindex, LoopVectorization.OperationStruct(0x00000000000000000000000000001536, 0x00000000000000000000000000000000, 0x00000000000000000000000000000000, 0x0000000000000000000000000008000a, 0x00000000000000000000000000000000, LoopVectorization.memload, 0x000b, 0x06), :LoopVectorization, :+, LoopVectorization.OperationStruct(0x00000000000000000000000000000025, 0x00000000000000000000000000000000, 0x00000000000000000000000000000000, 0x00000000000000000000000000020007, 0x00000000000000000000000000000000, LoopVectorization.compute, 0x000c, 0x00), :LoopVectorization, :+, LoopVectorization.OperationStruct(0x00000000000000000000000000000046, 0x00000000000000000000000000000000, 0x00000000000000000000000000000000, 0x00000000000000000000000000040009, 0x00000000000000000000000000000000, LoopVectorization.compute, 0x000d, 0x00), :LoopVectorization, :getindex, LoopVectorization.OperationStruct(0x00000000000000000000000000002546, 0x00000000000000000000000000000000, 0x00000000000000000000000000000000, 0x000000000000000000000000000c000d, 0x00000000000000000000000000000000, LoopVectorization.memload, 0x000e, 0x07), :LoopVectorization, :mul_fast, LoopVectorization.OperationStruct(0x00000000000000000000000000153624, 0x00000000000000000000000000000000, 0x00000000000000000000000000000000, 0x000000000000000000000000000b000e, 0x00000000000000000000000000000000, LoopVectorization.compute, 0x000f, 0x00), :LoopVectorization, :vfmadd_fast, LoopVectorization.OperationStruct(0x00000000000000000000000000561324, 0x00000000000000000000000000000056, 0x00000000000000000000000000000000, 0x000000000000000000000006000f0005, 0x00000000000000000000000000000000, LoopVectorization.compute, 0x0005, 0x00), :LoopVectorization, :identity, LoopVectorization.OperationStruct(0x00000000000000000000000000001234, 0x00000000000000000000000000000056, 0x00000000000000000000000000000000, 0x00000000000000000000000000000010, 0x00000000000000000000000000000000, LoopVectorization.compute, 0x0005, 0x00), :LoopVectorization, :getindex, LoopVectorization.OperationStruct(0x00000000000000000000000000000034, 0x00000000000000000000000000000000, 0x00000000000000000000000000000000, 0x00000000000000000000000000000000, 0x00000000000000000000000000000000, LoopVectorization.memload, 0x0010, 0x08), :numericconstant, Symbol("###reduction##zero###38###"), LoopVectorization.OperationStruct(0x00000000000000000000000000000034, 0x00000000000000000000000000000000, 0x00000000000000000000000000005612, 0x00000000000000000000000000000000, 0x00000000000000000000000000000000, LoopVectorization.constant, 0x0011, 0x00), :LoopVectorization, :add_fast, LoopVectorization.OperationStruct(0x00000000000000000000000000001234, 0x00000000000000000000000000005612, 0x00000000000000000000000000000000, 0x00000000000000000000000000110013, 0x00000000000000000000000000000000, LoopVectorization.compute, 0x0011, 0x00), :LoopVectorization, :reduced_add, LoopVectorization.OperationStruct(0x00000000000000000000000000000034, 0x00000000000000000000000000000012, 0x00000000000000000000000000000000, 0x00000000000000000000000000140012, 0x00000000000000000000000000000000, LoopVectorization.compute, 0x0010, 0x00), :LoopVectorization, :setindex!, LoopVectorization.OperationStruct(0x00000000000000000000000000000034, 0x00000000000000000000000000000000, 0x00000000000000000000000000000000, 0x00000000000000000000000000000015, 0x00000000000000000000000000000000, LoopVectorization.memstore, 0x0012, 0x08))}(), Val{(LoopVectorization.ArrayRefStruct{:range_a, Symbol("##vptr##_range_a")}(0x00000000000000000000000000000001, 0x00000000000000000000000000000001, 0x00000000000000000000000000000000, 0x00000000000000000000000000000001), LoopVectorization.ArrayRefStruct{:range_a, Symbol("##vptr##_range_a")}(0x00000000000000000000000000000001, 0x00000000000000000000000000000002, 0x00000000000000000000000000000000, 0x00000000000000000000000000000001), LoopVectorization.ArrayRefStruct{:range_b, Symbol("##vptr##_range_b")}(0x00000000000000000000000000000001, 0x00000000000000000000000000000003, 0x00000000000000000000000000000000, 0x00000000000000000000000000000001), LoopVectorization.ArrayRefStruct{:range_b, Symbol("##vptr##_range_b")}(0x00000000000000000000000000000001, 0x00000000000000000000000000000004, 0x00000000000000000000000000000000, 0x00000000000000000000000000000001), LoopVectorization.ArrayRefStruct{:src, Symbol("##vptr##_src")}(0x00000000000000000000000000000101, 0x00000000000000000000000000000506, 0x00000000000000000000000000000000, 0x00000000000000000000000000000101), LoopVectorization.ArrayRefStruct{:src, Symbol("##vptr##_src")}(0x00000000000000000000000000000202, 0x0000000000000000000000000000080a, 0x00000000000000000000000000000000, 0x00000000000000000000000000000101), LoopVectorization.ArrayRefStruct{:src, Symbol("##vptr##_src")}(0x00000000000000000000000000000202, 0x00000000000000000000000000000c0d, 0x00000000000000000000000000000000, 0x00000000000000000000000000000101), LoopVectorization.ArrayRefStruct{:target_arr, Symbol("##vptr##_target_arr")}(0x00000000000000000000000000000101, 0x00000000000000000000000000000304, 0x00000000000000000000000000000000, 0x00000000000000000000000000000101))}(), Val{(0, (), (), (), (), ((5, LoopVectorization.HardInt), (19, LoopVectorization.IntOrFloat)), ())}(), Val{(:a1i, :a2i, :b1i, :b2i, :i_a, :i_b)}(), Base.Val(Base.typeof(var"##vargsym#1186")), LoopVectorization.flatten_to_tuple(var"##vargsym#1186")...)
end, Symbol("#range_a#preserve#buffer#"), Symbol("#range_b#preserve#buffer#"), Symbol("#src#preserve#buffer#"), Symbol("#target_arr#preserve#buffer#")))
        nothing
    else
        begin
            #= logging.jl:350 =#
            let
                #= logging.jl:351 =#
                var"#2193#level" = Base.CoreLogging.Warn
                #= logging.jl:352 =#
                var"#2194#std_level" = Base.CoreLogging.convert(Base.CoreLogging.LogLevel, var"#2193#level")
                #= logging.jl:353 =#
                if var"#2194#std_level" >= Base.CoreLogging._min_enabled_level[]
                    #= logging.jl:354 =#
                    var"#2195#group" = :condense_loopset
                    #= logging.jl:355 =#
                    var"#2196#_module" = Main
                    #= logging.jl:356 =#
                    var"#2197#logger" = Base.CoreLogging.current_logger_for_env(var"#2194#std_level", var"#2195#group", var"#2196#_module")
                    #= logging.jl:357 =#
                    if !(var"#2197#logger" === Base.CoreLogging.nothing)
                        #= logging.jl:358 =#
                        var"#2198#id" = :Main_bae305b4
                        #= logging.jl:361 =#
                        if Base.CoreLogging._invoked_shouldlog(var"#2197#logger", var"#2193#level", var"#2196#_module", var"#2195#group", var"#2198#id")
                            #= logging.jl:362 =#
                            var"#2199#file" = "/home/graham/.julia/packages/LoopVectorization/KKbgb/src/condense_loopset.jl"
                            #= logging.jl:363 =#
                            var"#2200#line" = 825
                            #= logging.jl:364 =#
                            local var"#2202#msg", var"#2203#kwargs"
                            #= logging.jl:365 =#
                            begin
                                    #= logging.jl:326 =#
                                    let var"#2201#err" = nothing
                                        #= logging.jl:327 =#
                                        if var"#2201#err" === Base.CoreLogging.nothing
                                            #= logging.jl:328 =#
                                            var"#2202#msg" = "#= DIR/loopvec_test.jl:13 =#:\n`LoopVectorization.check_args` on your inputs failed; running fallback `@inbounds @fastmath` loop instead.\nUse `warn_check_args=false`, e.g. `@turbo warn_check_args=false ...`, to disable this warning."
                                            #= logging.jl:329 =#
                                            var"#2203#kwargs" = (; maxlog = 1)
                                            #= logging.jl:330 =#
                                            true
                                        else
                                            #= logging.jl:332 =#
                                            Base.CoreLogging.logging_error(var"#2197#logger", var"#2193#level", var"#2196#_module", var"#2195#group", var"#2198#id", var"#2199#file", var"#2200#line", var"#2201#err", false)
                                            #= logging.jl:333 =#
                                            false
                                        end
                                    end
                                end && Base.CoreLogging.handle_message(var"#2197#logger", var"#2193#level", var"#2202#msg", var"#2196#_module", var"#2195#group", var"#2198#id", var"#2199#file", var"#2200#line"; var"#2203#kwargs"...)
                        end
                    end
                end
                #= logging.jl:371 =#
                Base.CoreLogging.nothing
            end
        end
        begin
            $(Expr(:inbounds, true))
            local var"#2218#val" = for a1i = eachindex(range_a), a2i = eachindex(range_a), b1i = eachindex(range_b), b2i = eachindex(range_b)
                        #= DIR/loopvec_test.jl:17 =#
                        a1 = range_a[a1i]
                        #= DIR/loopvec_test.jl:18 =#
                        a2 = range_a[a2i]
                        #= DIR/loopvec_test.jl:19 =#
                        b1 = range_b[b1i]
                        #= DIR/loopvec_test.jl:19 =#
                        b2 = range_b[b2i]
                        #= DIR/loopvec_test.jl:20 =#
                        contribution = 0
                        #= DIR/loopvec_test.jl:21 =#
                        for i_a = padded_axis_a, i_b = padded_axis_b
                            #= DIR/loopvec_test.jl:22 =#
                            contribution = Base.FastMath.add_fast(contribution, Base.FastMath.mul_fast(src[i_a, i_b], src[Base.FastMath.add_fast(i_a, a1), Base.FastMath.add_fast(i_b, b1)], src[Base.FastMath.add_fast(i_a, a2), Base.FastMath.add_fast(i_b, b2)]))
                        end
                        #= DIR/loopvec_test.jl:24 =#
                        begin
                            #= fastmath.jl:116 =#
                            var"##1187" = target_arr
                            #= fastmath.jl:117 =#
                            (var"##1188", var"##1189") = (b1i, b2i)
                            #= fastmath.jl:118 =#
                            var"##1187"[var"##1188", var"##1189"] = Base.FastMath.add_fast(var"##1187"[var"##1188", var"##1189"], contribution)
                        end
                    end
            $(Expr(:inbounds, :pop))
            var"#2218#val"
        end
    end
end