
import sys

def count(inp):
	return len(inp.split())

print(count(sys.stdin.read()))