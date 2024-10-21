import json
import os
import subprocess
import glob
import time


def load_specification(spec_file):
    """Load the overall specification from a JSON file."""
    print(f"Loading specification from {spec_file}")
    with open(spec_file, "r") as f:
        spec = json.load(f)
    print(f"Loaded specification: {json.dumps(spec, indent=4)}")
    return spec


def generate_data(varied_param, value, fixed_params, output_base_dir):
    """Generate signature and sample data for a given parameter value using SLURM."""
    print(f"Generating data for {varied_param} = {value}")

    # Create the base directory for the varied parameter (e.g., data/param_being_varied)
    base_dir = os.path.join(output_base_dir, varied_param)
    os.makedirs(base_dir, exist_ok=True)
    print(f"Created base directory: {base_dir}")

    # Create the subfolder for this specific parameter value (e.g., data/param_being_varied/value)
    subfolder_name = f"{value}"
    subfolder_path = os.path.join(base_dir, subfolder_name)
    os.makedirs(subfolder_path, exist_ok=True)
    print(f"Created subfolder for parameter value: {subfolder_path}")

    # Update the parameter dictionary to include the varied parameter's value
    current_params = fixed_params.copy()
    current_params[varied_param] = value

    # Generate the signature file
    sig_file = os.path.join(subfolder_path, current_params["output_fasta"])
    sig_command = (
        f"./gen_sig {current_params['number_of_signatures']} "
        f"{current_params['signature_length']} {current_params['signature_length']} "
        f"{current_params['n_ratio_signatures']} > {sig_file}"
    )
    print(f"Signature command: {sig_command}")

    # Generate the sample file
    sample_file = os.path.join(subfolder_path, current_params["output_fastq"])
    sample_command = (
        f"./gen_sample {sig_file} {current_params['num_no_virus']} {current_params['num_with_virus']} "
        f"{current_params['min_viruses']} {current_params['max_viruses']} {current_params['sample_length']} "
        f"{current_params['sample_length']} {current_params['min_phred']} {current_params['max_phred']} "
        f"{current_params['n_ratio_samples']} > {sample_file}"
    )
    print(f"Sample command: {sample_command}")

    # Create a SLURM batch script
    slurm_script_path = os.path.join(subfolder_path, "job.slurm")
    with open(slurm_script_path, "w") as f:
        f.write("#!/bin/bash\n")
        f.write("#SBATCH --job-name=gen_data\n")
        f.write("#SBATCH --output=slurm-%j.out\n")
        f.write("#SBATCH --error=slurm-%j.err\n")
        f.write("#SBATCH --ntasks=1\n")
        f.write("#SBATCH --time=01:00:00\n")  # Adjust the time limit as needed
        f.write("\n")
        f.write(f"{sig_command}\n")
        f.write(f"{sample_command}\n")
    print(f"Created SLURM batch script: {slurm_script_path}")

    # Submit the SLURM job
    result = subprocess.run(
        f"sbatch {slurm_script_path}", shell=True, capture_output=True, text=True
    )
    if result.returncode != 0:
        print(f"Error submitting SLURM job: {result.stderr}")
        return None

    # Save the SLURM job ID for later tracking
    job_id = result.stdout.strip().split()[-1]
    print(f"Submitted SLURM job. Job ID: {job_id}")
    with open(os.path.join(subfolder_path, "slurm_job_id.txt"), "w") as f:
        f.write(job_id)

    # Generate the subfolder's specification.json
    with open(os.path.join(subfolder_path, "specification.json"), "w") as f:
        json.dump(current_params, f, indent=4)
    print(
        f"Saved current specification to {os.path.join(subfolder_path, 'specification.json')}"
    )

    # Log the generation process
    with open(os.path.join(subfolder_path, "generation.log"), "w") as f:
        f.write(f"Generated data with {varied_param} = {value}\n")
        f.write(f"Signature command: {sig_command}\n")
        f.write(f"Sample command: {sample_command}\n")
        f.write(f"SLURM job ID: {job_id}\n")

    print(
        f"Data generation completed for {varied_param} = {value}. Logged details in {os.path.join(subfolder_path, 'generation.log')}"
    )
    return job_id


def wait_for_jobs_to_complete(job_ids):
    """Wait for all SLURM jobs to complete."""
    print(f"Waiting for {len(job_ids)} jobs to complete: {job_ids}")
    while job_ids:
        # Check the status of each job
        for job_id in job_ids[:]:  # Create a copy of the list to iterate over
            result = subprocess.run(
                f"squeue -j {job_id}", shell=True, capture_output=True, text=True
            )
            if result.returncode == 0 and job_id not in result.stdout:
                # Job is not in the queue, it has finished
                job_ids.remove(job_id)
                print(f"Job {job_id} has completed.")
            else:
                print(f"Job {job_id} is still running.")

        # Wait for a few seconds before checking again
        time.sleep(1)
    print("All jobs have completed.")


def process_specification_file(spec_file, output_dir):
    """Process a single specification file."""
    print(f"Processing specification file: {spec_file}")
    # Load the specification
    spec = load_specification(spec_file)

    # Identify the varied parameter
    varied_params = {
        k: v for k, v in spec.items() if isinstance(v, dict) and "range" in v
    }

    # Validate that there is only one varied parameter
    if len(varied_params) != 1:
        raise ValueError(
            f"There should be exactly one varied parameter in the specification {spec_file}."
        )
    varied_param, param_config = next(iter(varied_params.items()))
    print(f"Identified varied parameter: {varied_param}")

    # Extract the range and step for the varied parameter
    param_range = param_config["range"]
    param_step = param_config["step"]
    print(f"Parameter range: {param_range}, step: {param_step}")

    # Get the fixed parameters
    fixed_params = {k: v for k, v in spec.items() if k not in varied_params}
    print(f"Fixed parameters: {json.dumps(fixed_params, indent=4)}")

    # Loop over the specified range for the varied parameter
    start, end = param_range
    current_value = start
    job_ids = []
    while current_value <= end:
        # Generate data for the current value of the varied parameter
        print(f"Generating data for {varied_param} with value {current_value}")
        job_id = generate_data(varied_param, current_value, fixed_params, output_dir)
        if job_id:
            job_ids.append(job_id)
        current_value += param_step

    # Save the overall specification in the output directory
    spec_output_dir = os.path.join(output_dir)
    os.makedirs(spec_output_dir, exist_ok=True)
    with open(os.path.join(spec_output_dir, "overall_specification.json"), "w") as f:
        json.dump(spec, f, indent=4)
    print(
        f"Saved overall specification to {os.path.join(spec_output_dir, 'overall_specification.json')}"
    )

    # Return the list of job IDs for this specification
    return job_ids


def main():
    # Directory where specification files are stored
    spec_dir = "specifications"
    # Output directory for generated data
    output_dir = "data"

    # Get a list of all JSON specification files in the specifications folder
    spec_files = glob.glob(os.path.join(spec_dir, "*.json"))

    # List to collect all job IDs
    all_job_ids = []

    # Process each specification file
    for spec_file in spec_files:
        print(f"Starting processing for specification: {spec_file}")
        job_ids = process_specification_file(spec_file, output_dir)
        all_job_ids.extend(job_ids)
        print(f"Finished processing for specification: {spec_file}")

    # Wait for all jobs to complete
    wait_for_jobs_to_complete(all_job_ids)
    print("All jobs for all specifications have been completed.")


if __name__ == "__main__":
    main()
