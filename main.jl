#!/usr/bin/env julia

using Nas2Step

function main(nas_file; is_verify=false)
    println("Input NAS file: $nas_file")

    step_path = replace(nas_file, ".nas" => "_auto.step")
    step_file = Nas2Step.nas_to_step(nas_file, step_path=step_path)

     # Verify the output STEP file
    if is_verify
        println("\nVerifying output STEP file:")
        Nas2Step.verify_step(step_file)
    end

end

# Run the example with a realistic volume-only NAS file from input Arguments
if length(ARGS) > 0
    nas_file = ARGS[1]
    main(nas_file; is_verify=false)
else
    println("No input NAS file specified.")
end