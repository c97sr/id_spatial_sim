# this has been ported to ebola_build.exe
# leaving here for future reference
function randint(n)
{
    return int(n * rand())
}

BEGIN {
    FS=","
    OFS=","
    i=0
    max=0
}

# this runs on the first input file to awk
FNR==NR {
    demographic[i]=$2","$3","$4
    number[i]=$5
    i+=1
    max+=$5
    next
}

# this runs on the second input file
{
    if ($5==hhs) {
        j=0
        sum=0
        rnd=randint(max)+1
        while (sum < rnd) {
            sum+=number[j]
            j++
        }
        print $1,$2,$3,demographic[j-1]
    }
}

