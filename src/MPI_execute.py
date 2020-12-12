import subprocess

cmd = "main_mpi.cpp"
# Example 
# cmd = HelloWorld.c 
for i in range(0,100):
    subprocess.call(["mpicxx", '-std=c++14', cmd, '-o', "main"])  # For Compiling
    subprocess.call(["mpiexec",'--oversubscribe', "-np", '4', "./main"])