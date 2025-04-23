import subprocess
import sys

def execute_command(cmd, shell=False):
    if shell:
        process = subprocess.Popen(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            bufsize=1,
            shell=True
        )
    # shell=False: the command is executed directly without going through the system's command shell (like bash or cmd.exe). 
    # This is generally safer and avoids potential security issues related to shell injection. 
    # The command must be passed as a list in this case.
    elif not shell:
        process = subprocess.Popen(
            cmd,
            # both normal output and error messages will be captured together through process.stdout
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            bufsize=1
        )
    else:
        raise ValueError("incorrect shell argument")
    
    # Read the output line by line in real-time
    for line in iter(process.stdout.readline, ''):
        print(line, end='')  # Print without extra newline
        sys.stdout.flush()   # Make sure it's displayed immediately
    # Wait for the process to complete
    process.communicate()
    return process.returncode