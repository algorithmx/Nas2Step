
# ============================================================================
# Helper functions for concise reporting
# ============================================================================

function fmt_pct(n, total)
    return round(100 * n / max(total, 1), digits=1)
end

function fmt_stat(x; digits=4)
    return round(x, sigdigits=digits)
end

function print_volume_stats(vs)
    println("  Volume: [$(fmt_stat(vs.min)), $(fmt_stat(vs.max))], " *
           "mean=$(fmt_stat(vs.mean)), std=$(fmt_stat(vs.std))")
end

function print_inverted_summary(vol_check)
    inv_pct = fmt_pct(vol_check.inverted_count, vol_check.total_elements)
    println("  $(vol_check.inverted_count)/$(vol_check.total_elements) inverted ($(inv_pct)%), " *
           "$(vol_check.zero_volume_count) zero-volume")
end

function print_inverted_details(inverted_elements, max_show=10)
    n_show = min(max_show, length(inverted_elements))
    sorted = sort(inverted_elements, by = e -> e.volume)
    
    for i in 1:n_show
        elem = sorted[i]
        println("  [$(i)] Elem $(elem.element_id): V=$(fmt_stat(elem.volume, digits=5)), " *
               "nodes=$(elem.node_ids), center=$(round.(elem.centroid, digits=2))")
    end
    
    if length(inverted_elements) > n_show
        println("  ... $(length(inverted_elements) - n_show) more")
    end
end

function print_swap_test_result(swap_test)
    println("  🔬 Swap test (1⇄2): $(swap_test.original_inverted) → $(swap_test.swapped_inverted) " *
           "[$(swap_test.fixed_count) fixed, $(fmt_pct(swap_test.improvement_ratio, 1))%]")
    if swap_test.would_fix
        println("  ✅ Would fix! This confirms convention mismatch.")
    else
        println("  ⚠️  Wouldn't fully fix. May have real defects.")
    end
end



# Simple JSON writer (no external deps)
function write_json(io::IO, obj, indent::Int)
    ind = "  " ^ indent
    if obj isa Dict
        println(io, "{")
        keys_list = collect(keys(obj))
        for (i, k) in enumerate(keys_list)
            print(io, ind, "  \"", k, "\": ")
            write_json(io, obj[k], indent + 1)
            println(io, i < length(keys_list) ? "," : "")
        end
        print(io, ind, "}")
    elseif obj isa AbstractVector
        println(io, "[")
        for (i, v) in enumerate(obj)
            print(io, ind, "  ")
            write_json(io, v, indent + 1)
            println(io, i < length(obj) ? "," : "")
        end
        print(io, ind, "]")
    elseif obj isa AbstractString
        print(io, "\"", obj, "\"")
    elseif obj isa Number
        print(io, obj)
    elseif obj isa Bool
        print(io, obj ? "true" : "false")
    elseif obj === nothing
        print(io, "null")
    else
        print(io, "\"", string(obj), "\"")
    end
end
