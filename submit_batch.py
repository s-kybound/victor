import os
import subprocess
import time
import re
import csv
from concurrent.futures import ThreadPoolExecutor


def create_slurm_script(
    folder_path, gpu_type, script_name, matcher_command, log_filename
):
    """Create a SLURM batch script for the specified GPU type."""
    # Set GPU type and constraints based on whether it's H100 or A100
    if gpu_type == "h100":
        gpu_config = "h100-96"
        constraint = "xgpi"
    elif gpu_type == "a100":
        gpu_config = "a100-40"
        constraint = "xgph"
    else:
        raise ValueError(f"Unknown GPU type: {gpu_type}")

    slurm_template = f"""#!/bin/bash
#SBATCH --job-name={gpu_type}_matcher_job          # Job name
#SBATCH --output={log_filename}                    # Output file for logging
#SBATCH --ntasks=1                                 # Number of tasks
#SBATCH --cpus-per-task=1                          # Number of CPU cores per task
#SBATCH --mem=20G                                  # Total memory
#SBATCH --gpus={gpu_config}                        # GPU type and number
#SBATCH --constraint={constraint}                  # GPU constraint
#SBATCH --time=00:10:00                            # Time limit

# Run the matcher command directly
{matcher_command}
"""
    script_path = os.path.join(folder_path, script_name)
    with open(script_path, "w") as script_file:
        script_file.write(slurm_template)
    return script_path


def submit_job(script_path):
    """Submit a batch job and return the job ID."""
    print(f"Submitting job with script: {script_path}")
    result = subprocess.run(["sbatch", script_path], capture_output=True, text=True)

    # Parse the job ID from the sbatch output
    if result.returncode != 0:
        print(f"Failed to submit job: {result.stderr}")
        return None

    job_id = re.search(r"Submitted batch job (\d+)", result.stdout)
    if not job_id:
        print("Could not parse job ID from sbatch output.")
        return None

    job_id = job_id.group(1)
    print(f"Submitted job with ID: {job_id}")
    return job_id


def monitor_jobs_async(job_ids):
    """Monitor the status of jobs asynchronously."""
    pending_jobs = set(job_ids)
    completed_jobs = set()

    print(f"Monitoring {len(job_ids)} jobs asynchronously...")
    while pending_jobs:
        for job_id in list(pending_jobs):
            # Check job status using squeue
            check = subprocess.run(
                ["squeue", "--job", job_id], capture_output=True, text=True
            )
            if job_id not in check.stdout:
                print(f"Job {job_id} has completed.")
                pending_jobs.remove(job_id)
                completed_jobs.add(job_id)
        time.sleep(1)  # Check every second

    print("All jobs have completed.")
    return completed_jobs


def extract_timing_from_log(log_file):
    """Extract the timing from the log file."""
    try:
        with open(log_file, "r") as f:
            log_content = f.read()
        match = re.search(
            r"\(FOR AUTOMATED CHECKING\) Total runMatcher time:(\d+\.\d+)", log_content
        )
        if match:
            timing = float(match.group(1))
            print(f"Extracted timing from {log_file}: {timing}s")
            return timing
        else:
            print(f"Timing not found in {log_file}")
    except FileNotFoundError:
        print(f"Log file not found: {log_file}")
    return None


def run_batch_jobs_async(folder_path, k, matcher_command_h100, matcher_command_a100):
    """Run batch jobs k times asynchronously and log the timings."""
    h100_jobs = []
    a100_jobs = []
    h100_timings = []
    a100_timings = []

    # Submit all jobs asynchronously
    with ThreadPoolExecutor(max_workers=4) as executor:
        for i in range(k):
            # Submit H100 job
            h100_log = os.path.join(folder_path, f"h100_output_{i + 1}.txt")
            h100_script = create_slurm_script(
                folder_path,
                "h100",
                f"h100_matcher_job_{i + 1}.sh",
                matcher_command_h100,
                h100_log,
            )
            h100_job_id = executor.submit(submit_job, h100_script).result()
            if h100_job_id:
                h100_jobs.append((h100_job_id, h100_log))

            # Submit A100 job
            a100_log = os.path.join(folder_path, f"a100_output_{i + 1}.txt")
            a100_script = create_slurm_script(
                folder_path,
                "a100",
                f"a100_matcher_job_{i + 1}.sh",
                matcher_command_a100,
                a100_log,
            )
            a100_job_id = executor.submit(submit_job, a100_script).result()
            if a100_job_id:
                a100_jobs.append((a100_job_id, a100_log))

    # Monitor all jobs asynchronously
    h100_job_ids = [job_id for job_id, _ in h100_jobs]
    a100_job_ids = [job_id for job_id, _ in a100_jobs]
    monitor_jobs_async(h100_job_ids + a100_job_ids)

    # Extract timings once all jobs have completed
    for _, h100_log in h100_jobs:
        h100_time = extract_timing_from_log(h100_log)
        if h100_time is not None:
            h100_timings.append(h100_time)

    for _, a100_log in a100_jobs:
        a100_time = extract_timing_from_log(a100_log)
        if a100_time is not None:
            a100_timings.append(a100_time)

    # Return average timings
    avg_h100_timing = sum(h100_timings) / len(h100_timings) if h100_timings else None
    avg_a100_timing = sum(a100_timings) / len(a100_timings) if a100_timings else None

    return avg_h100_timing, avg_a100_timing


def main(output_dir, k):
    """Main function to handle batch job submissions and timing collection."""
    results = [("param_being_varied", "param_value", "h100_timing", "a100_timing")]

    # Traverse the folder structure data/<param_being_varied>/<value>
    for param_being_varied in os.listdir(output_dir):
        param_dir = os.path.join(output_dir, param_being_varied)
        if not os.path.isdir(param_dir):
            continue

        for folder_name in os.listdir(param_dir):
            folder_path = os.path.join(param_dir, folder_name)
            if not os.path.isdir(folder_path):
                continue

            # Extract the parameter value
            param_value = folder_name

            # Paths to samp.fastq and sig.fasta
            samp_path = os.path.join(folder_path, "samp.fastq")
            sig_path = os.path.join(folder_path, "sig.fasta")

            # Matcher commands for H100 and A100
            matcher_command_h100 = f"./matcher-h100 {samp_path} {sig_path}"
            matcher_command_a100 = f"./matcher-a100 {samp_path} {sig_path}"

            # Run batch jobs asynchronously and collect average timings
            print(f"Running batch jobs for {param_being_varied} = {param_value}")
            avg_h100_timing, avg_a100_timing = run_batch_jobs_async(
                folder_path, k, matcher_command_h100, matcher_command_a100
            )

            # Append the results
            results.append(
                (param_being_varied, param_value, avg_h100_timing, avg_a100_timing)
            )

    # Write the results to a CSV file
    csv_path = os.path.join(output_dir, "timings.csv")
    print(f"Writing results to {csv_path}")
    with open(csv_path, "w", newline="") as csvfile:
        csvwriter = csv.writer(csvfile)
        csvwriter.writerows(results)


if __name__ == "__main__":
    output_dir = "data"  # Change this to your actual output directory path
    k = 3  # Number of repetitions for each job
    main(output_dir, k)
