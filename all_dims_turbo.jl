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
        #= DIR/loopvec_test.jl:13 =#
        nothing
    end
    var"###looprangeb1i###1###" = LoopVectorization.canonicalize_range(begin
                $(Expr(:inbounds, true))
                local var"#2219#val" = eachindex(range_b)
                $(Expr(:inbounds, :pop))
                var"#2219#val"
            end)
    var"###looplenb1i###2###" = ArrayInterface.static_length(var"###looprangeb1i###1###")
    var"###b1i_loop_lower_bound###3###" = LoopVectorization.maybestaticfirst(var"###looprangeb1i###1###")
    var"###b1i_loop_upper_bound###4###" = LoopVectorization.maybestaticlast(var"###looprangeb1i###1###")
    var"###looprangeb2i###5###" = LoopVectorization.canonicalize_range(begin
                $(Expr(:inbounds, true))
                local var"#2220#val" = eachindex(range_b)
                $(Expr(:inbounds, :pop))
                var"#2220#val"
            end)
    var"###looplenb2i###6###" = ArrayInterface.static_length(var"###looprangeb2i###5###")
    var"###b2i_loop_lower_bound###7###" = LoopVectorization.maybestaticfirst(var"###looprangeb2i###5###")
    var"###b2i_loop_upper_bound###8###" = LoopVectorization.maybestaticlast(var"###looprangeb2i###5###")
    var"###looprangei_a###10###" = LoopVectorization.canonicalize_range(begin
                $(Expr(:inbounds, true))
                local var"#2221#val" = padded_axis_a
                $(Expr(:inbounds, :pop))
                var"#2221#val"
            end)
    var"###loopleni_a###11###" = ArrayInterface.static_length(var"###looprangei_a###10###")
    var"###i_a_loop_lower_bound###12###" = LoopVectorization.maybestaticfirst(var"###looprangei_a###10###")
    var"###i_a_loop_upper_bound###13###" = LoopVectorization.maybestaticlast(var"###looprangei_a###10###")
    var"###i_a_loop_step###14###" = LoopVectorization.static_step(var"###looprangei_a###10###")
    var"###looprangei_b###15###" = LoopVectorization.canonicalize_range(begin
                $(Expr(:inbounds, true))
                local var"#2222#val" = padded_axis_b
                $(Expr(:inbounds, :pop))
                var"#2222#val"
            end)
    var"###loopleni_b###16###" = ArrayInterface.static_length(var"###looprangei_b###15###")
    var"###i_b_loop_lower_bound###17###" = LoopVectorization.maybestaticfirst(var"###looprangei_b###15###")
    var"###i_b_loop_upper_bound###18###" = LoopVectorization.maybestaticlast(var"###looprangei_b###15###")
    var"###i_b_loop_step###19###" = LoopVectorization.static_step(var"###looprangei_b###15###")
    if LoopVectorization.check_args(range_b, src, target_arr)
        (var"##vptr##_range_b", var"#range_b#preserve#buffer#") = LoopVectorization.stridedpointer_preserve(range_b)
        (var"##vptr##_src", var"#src#preserve#buffer#") = LoopVectorization.stridedpointer_preserve(src)
        (var"##vptr##_target_arr", var"#target_arr#preserve#buffer#") = LoopVectorization.stridedpointer_preserve(target_arr)
        var"####grouped#strided#pointer####27###" = Core.getfield(LoopVectorization.grouped_strided_pointer((LoopVectorization.densewrapper(LoopVectorization.gespf1(var"##vptr##_range_b", (LoopVectorization.similardims(var"###looprangeb1i###1###", VectorizationBase.NullStep()),)), range_b), LoopVectorization.densewrapper(LoopVectorization.gespf1(var"##vptr##_src", (var"###i_a_loop_lower_bound###12###", LoopVectorization.similardims(var"###looprangei_b###15###", VectorizationBase.NullStep()))), src), LoopVectorization.densewrapper(LoopVectorization.gespf1(var"##vptr##_target_arr", (LoopVectorization.similardims(var"###looprangeb1i###1###", VectorizationBase.NullStep()), LoopVectorization.similardims(var"###looprangeb2i###5###", VectorizationBase.NullStep()))), target_arr)), Val{()}()), 1)
        $(Expr(:gc_preserve, quote
    var"##vargsym#1190" = ((var"###looprangeb1i###1###", var"###looprangeb2i###5###", LoopVectorization.zerorangestart(var"###looprangei_a###10###"), var"###looprangei_b###15###"), (var"####grouped#strided#pointer####27###",))
    var"##Tloopeltype##" = LoopVectorization.promote_type(LoopVectorization.eltype(range_b), LoopVectorization.eltype(src), LoopVectorization.eltype(target_arr))
    var"##Wvecwidth##" = LoopVectorization.pick_vector_width(var"##Tloopeltype##")
    LoopVectorization._turbo_!(LoopVectorization.avx_config_val(Val{(false, 0, 0, 0, false, 0x0000000000000001)}(), var"##Wvecwidth##"), Val{(:LoopVectorization, :getindex, LoopVectorization.OperationStruct(0x00000000000000000000000000000001, 0x00000000000000000000000000000000, 0x00000000000000000000000000000000, 0x00000000000000000000000000000000, 0x00000000000000000000000000000000, LoopVectorization.memload, 0x0001, 0x01), :LoopVectorization, :getindex, LoopVectorization.OperationStruct(0x00000000000000000000000000000002, 0x00000000000000000000000000000000, 0x00000000000000000000000000000000, 0x00000000000000000000000000000000, 0x00000000000000000000000000000000, LoopVectorization.memload, 0x0002, 0x02), :numericconstant, Symbol("###0###9###"), LoopVectorization.OperationStruct(0x00000000000000000000000000000012, 0x00000000000000000000000000000000, 0x00000000000000000000000000000034, 0x00000000000000000000000000000000, 0x00000000000000000000000000000000, LoopVectorization.constant, 0x0003, 0x00), :LoopVectorization, :getindex, LoopVectorization.OperationStruct(0x00000000000000000000000000000034, 0x00000000000000000000000000000000, 0x00000000000000000000000000000000, 0x00000000000000000000000000000000, 0x00000000000000000000000000000000, LoopVectorization.memload, 0x0004, 0x03), :i_b, :i_b, LoopVectorization.OperationStruct(0x00000000000000000000000000000004, 0x00000000000000000000000000000000, 0x00000000000000000000000000000000, 0x00000000000000000000000000000000, 0x00000000000000000000000000000000, LoopVectorization.loopvalue, 0x0005, 0x00), :LoopVectorization, :+, LoopVectorization.OperationStruct(0x00000000000000000000000000000014, 0x00000000000000000000000000000000, 0x00000000000000000000000000000000, 0x00000000000000000000000000010005, 0x00000000000000000000000000000000, LoopVectorization.compute, 0x0006, 0x00), :LoopVectorization, :getindex, LoopVectorization.OperationStruct(0x00000000000000000000000000000314, 0x00000000000000000000000000000000, 0x00000000000000000000000000000000, 0x00000000000000000000000000000006, 0x00000000000000000000000000000000, LoopVectorization.memload, 0x0007, 0x04), :LoopVectorization, :+, LoopVectorization.OperationStruct(0x00000000000000000000000000000024, 0x00000000000000000000000000000000, 0x00000000000000000000000000000000, 0x00000000000000000000000000020005, 0x00000000000000000000000000000000, LoopVectorization.compute, 0x0008, 0x00), :LoopVectorization, :getindex, LoopVectorization.OperationStruct(0x00000000000000000000000000000324, 0x00000000000000000000000000000000, 0x00000000000000000000000000000000, 0x00000000000000000000000000000008, 0x00000000000000000000000000000000, LoopVectorization.memload, 0x0009, 0x05), :LoopVectorization, :mul_fast, LoopVectorization.OperationStruct(0x00000000000000000000000000003142, 0x00000000000000000000000000000000, 0x00000000000000000000000000000000, 0x00000000000000000000000000070009, 0x00000000000000000000000000000000, LoopVectorization.compute, 0x000a, 0x00), :LoopVectorization, :vfmadd_fast, LoopVectorization.OperationStruct(0x00000000000000000000000000003412, 0x00000000000000000000000000000034, 0x00000000000000000000000000000000, 0x000000000000000000000004000a0003, 0x00000000000000000000000000000000, LoopVectorization.compute, 0x0003, 0x00), :LoopVectorization, :identity, LoopVectorization.OperationStruct(0x00000000000000000000000000000012, 0x00000000000000000000000000000034, 0x00000000000000000000000000000000, 0x0000000000000000000000000000000b, 0x00000000000000000000000000000000, LoopVectorization.compute, 0x0003, 0x00), :LoopVectorization, :getindex, LoopVectorization.OperationStruct(0x00000000000000000000000000000012, 0x00000000000000000000000000000000, 0x00000000000000000000000000000000, 0x00000000000000000000000000000000, 0x00000000000000000000000000000000, LoopVectorization.memload, 0x000b, 0x06), :LoopVectorization, :add_fast, LoopVectorization.OperationStruct(0x00000000000000000000000000000012, 0x00000000000000000000000000000034, 0x00000000000000000000000000000000, 0x000000000000000000000000000c000d, 0x00000000000000000000000000000000, LoopVectorization.compute, 0x000b, 0x00), :LoopVectorization, :setindex!, LoopVectorization.OperationStruct(0x00000000000000000000000000000012, 0x00000000000000000000000000000034, 0x00000000000000000000000000000000, 0x0000000000000000000000000000000e, 0x00000000000000000000000000000000, LoopVectorization.memstore, 0x000c, 0x06))}(), Val{(LoopVectorization.ArrayRefStruct{:range_b, Symbol("##vptr##_range_b")}(0x00000000000000000000000000000001, 0x00000000000000000000000000000001, 0x00000000000000000000000000000000, 0x00000000000000000000000000000001), LoopVectorization.ArrayRefStruct{:range_b, Symbol("##vptr##_range_b")}(0x00000000000000000000000000000001, 0x00000000000000000000000000000002, 0x00000000000000000000000000000000, 0x00000000000000000000000000000001), LoopVectorization.ArrayRefStruct{:src, Symbol("##vptr##_src")}(0x00000000000000000000000000000101, 0x00000000000000000000000000000304, 0x00000000000000000000000000000000, 0x00000000000000000000000000000101), LoopVectorization.ArrayRefStruct{:src, Symbol("##vptr##_src")}(0x00000000000000000000000000000102, 0x00000000000000000000000000000306, 0x00000000000000000000000000000000, 0x00000000000000000000000000000101), LoopVectorization.ArrayRefStruct{:src, Symbol("##vptr##_src")}(0x00000000000000000000000000000102, 0x00000000000000000000000000000308, 0x00000000000000000000000000000000, 0x00000000000000000000000000000101), LoopVectorization.ArrayRefStruct{:target_arr, Symbol("##vptr##_target_arr")}(0x00000000000000000000000000000101, 0x00000000000000000000000000000102, 0x00000000000000000000000000000000, 0x00000000000000000000000000000101))}(), Val{(0, (), (), (), (), ((3, LoopVectorization.HardInt),), ())}(), Val{(:b1i, :b2i, :i_a, :i_b)}(), Base.Val(Base.typeof(var"##vargsym#1190")), LoopVectorization.flatten_to_tuple(var"##vargsym#1190")...)
end, Symbol("#range_b#preserve#buffer#"), Symbol("#src#preserve#buffer#"), Symbol("#target_arr#preserve#buffer#")))
        nothing
    else
        begin
            #= logging.jl:350 =#
            let
                #= logging.jl:351 =#
                var"#2223#level" = Base.CoreLogging.Warn
                #= logging.jl:352 =#
                var"#2224#std_level" = Base.CoreLogging.convert(Base.CoreLogging.LogLevel, var"#2223#level")
                #= logging.jl:353 =#
                if var"#2224#std_level" >= Base.CoreLogging._min_enabled_level[]
                    #= logging.jl:354 =#
                    var"#2225#group" = :condense_loopset
                    #= logging.jl:355 =#
                    var"#2226#_module" = Main
                    #= logging.jl:356 =#
                    var"#2227#logger" = Base.CoreLogging.current_logger_for_env(var"#2224#std_level", var"#2225#group", var"#2226#_module")
                    #= logging.jl:357 =#
                    if !(var"#2227#logger" === Base.CoreLogging.nothing)
                        #= logging.jl:358 =#
                        var"#2228#id" = :Main_bae305b5
                        #= logging.jl:361 =#
                        if Base.CoreLogging._invoked_shouldlog(var"#2227#logger", var"#2223#level", var"#2226#_module", var"#2225#group", var"#2228#id")
                            #= logging.jl:362 =#
                            var"#2229#file" = "/home/graham/.julia/packages/LoopVectorization/KKbgb/src/condense_loopset.jl"
                            #= logging.jl:363 =#
                            var"#2230#line" = 825
                            #= logging.jl:364 =#
                            local var"#2232#msg", var"#2233#kwargs"
                            #= logging.jl:365 =#
                            begin
                                    #= logging.jl:326 =#
                                    let var"#2231#err" = nothing
                                        #= logging.jl:327 =#
                                        if var"#2231#err" === Base.CoreLogging.nothing
                                            #= logging.jl:328 =#
                                            var"#2232#msg" = "#= DIR/loopvec_test.jl:13 =#:\n`LoopVectorization.check_args` on your inputs failed; running fallback `@inbounds @fastmath` loop instead.\nUse `warn_check_args=false`, e.g. `@turbo warn_check_args=false ...`, to disable this warning."
                                            #= logging.jl:329 =#
                                            var"#2233#kwargs" = (; maxlog = 1)
                                            #= logging.jl:330 =#
                                            true
                                        else
                                            #= logging.jl:332 =#
                                            Base.CoreLogging.logging_error(var"#2227#logger", var"#2223#level", var"#2226#_module", var"#2225#group", var"#2228#id", var"#2229#file", var"#2230#line", var"#2231#err", false)
                                            #= logging.jl:333 =#
                                            false
                                        end
                                    end
                                end && Base.CoreLogging.handle_message(var"#2227#logger", var"#2223#level", var"#2232#msg", var"#2226#_module", var"#2225#group", var"#2228#id", var"#2229#file", var"#2230#line"; var"#2233#kwargs"...)
                        end
                    end
                end
                #= logging.jl:371 =#
                Base.CoreLogging.nothing
            end
        end
        begin
            $(Expr(:inbounds, true))
            local var"#2248#val" = for b1i = eachindex(range_b), b2i = eachindex(range_b)
                        #= DIR/loopvec_test.jl:19 =#
                        b1 = range_b[b1i]
                        #= DIR/loopvec_test.jl:19 =#
                        b2 = range_b[b2i]
                        #= DIR/loopvec_test.jl:20 =#
                        contribution = 0
                        #= DIR/loopvec_test.jl:21 =#
                        for i_a = padded_axis_a, i_b = padded_axis_b
                            #= DIR/loopvec_test.jl:22 =#
                            contribution = Base.FastMath.add_fast(contribution, Base.FastMath.mul_fast(src[i_a, i_b], src[i_a, Base.FastMath.add_fast(i_b, b1)], src[i_a, Base.FastMath.add_fast(i_b, b2)]))
                        end
                        #= DIR/loopvec_test.jl:24 =#
                        begin
                            #= fastmath.jl:116 =#
                            var"##1191" = target_arr
                            #= fastmath.jl:117 =#
                            (var"##1192", var"##1193") = (b1i, b2i)
                            #= fastmath.jl:118 =#
                            var"##1191"[var"##1192", var"##1193"] = Base.FastMath.add_fast(var"##1191"[var"##1192", var"##1193"], contribution)
                        end
                    end
            $(Expr(:inbounds, :pop))
            var"#2248#val"
        end
    end
end