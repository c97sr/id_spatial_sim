
function randint(n)
{
    return int(n * rand())
}

function uniform(min,max)
{
    return min+randint(max-min)
}

BEGIN {
    FS=","
    OFS=","
}

{
    for(i=0;i<$4;i++) print $1,uniform(0,19)
    for(i=0;i<$5;i++) print $1,uniform(20,64)
    for(i=0;i<$6;i++) print $1,uniform(65,85)
}

