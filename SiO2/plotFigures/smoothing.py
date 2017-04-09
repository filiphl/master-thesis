from numpy import zeros

def smoother(degree, function):
    smoothed = zeros(len(function))
    smoothed[0] = function[0]

    for i in xrange(1,degree+1):
        smoothed[i] = smooth(i,function, i)
    for i in xrange(degree, len(function)-degree):
        smoothed[i] = smooth(degree,function,i)
    for i in xrange(2,degree+2):
        smoothed[-i] = smooth(i-1,function, -i)
    smoothed[-1] = function[-1]
    return smoothed

def smooth(degree, function, atIndex):
    value            = 0.0
    dividor          = 0.0
    localCoeffisient = 1.0
    for i in xrange(-degree, degree+1):
        dividor += localCoeffisient
        value += localCoeffisient*function[atIndex+i]
        if i < 0:
            localCoeffisient += 1
        if i == 0:
            localCoeffisient -= 1
        if i > 0:
            localCoeffisient -= 1

    localCoeffisient +=1
    return value/dividor
