import subprocess

a = subprocess.run(["M2", "--script", "test.m2", "8"],capture_output=True).stdout
# os.system("M2 --script test.m2 8")
print((str(a)))