
#
with open("../../Data/raw_counts.txt", 'r') as w:
    w = w.readlines()
    count = 0
    for i in w[0].split():
        count += 1
    print ("RAW counts patient ID: {}".format(count-1))

    
with open("../data/raw_counts_NASH_code.txt", 'r') as w:
    w = w.readlines()
    count = 0
    z = w[0]
    # print (z.split())
    for i in w[0].split():
        count += 1
    print ("RAW counts NASH code: {}".format(count-1))
