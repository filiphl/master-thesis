import re
from sys import argv

infile  = open("system.prepare", "r")
content = infile.read()

if len(argv)<2:
    print "Use a unit cell length as input."
    exit()
else:
    a = float(argv[1])

'''
sphere = re.findall(r'sphere\s(\d+\.\d+)\s(\d+\.\d+)\s(\d+\.\d+)\s(\d+\.\d+)',content)
sphere = [float(match) for match in sphere[0]]
print sphere

myre = re.escape(str(a))
content2 = re.sub(r'sphere\s(\d+\.\d+)\s(\d+\.\d+)\s(\d+\.\d+)\s(\d+\.\d+)', (str(a))(3), content)
'''

numbers = re.findall(r'(\d+\.\d+)',content)
numbers = [float(match) for match in numbers]


def twelveToOneSixSix(match):
    match = float(match.group())
    return str(match*7.166/7.12)

def oneSixSixToTwelve(match):
    match = float(match.group())
    return str(match*7.12/7.166)


if a == 7.12:
    content2 = re.sub(r'(\d+\.\d+)', oneSixSixToTwelve, content)

if a == 7.166:
    content2 = re.sub(r'(\d+\.\d+)', twelveToOneSixSix, content)

print content2
