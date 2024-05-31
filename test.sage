import subprocess

# a = subprocess.run(["M2", "--script", "test.m2", "81","a9"],capture_output=True).stdout
# os.system("M2 --script test.m2 8")
# print((str(a)))

def run_macualay2_program(function_name, args):
	file_name = "compute_griffiths.m2"
	ls = ["M2", "--script", file_name, function_name]
	for arg in args:
		ls.append(str(arg))
	a = subprocess.run(ls,capture_output=True).stdout
	return a

a = run_macualay2_program("computeGriffithsRing", [2,"[x_0^3+x_1^3 + x_2^3]"])

print(str(a))