set -x

function show()
{
  f=$1
  en1=`sfget parform=n n1 < $f`
  en2=`sfget parform=n n2 < $f`
  en3=`sfget parform=n n3 < $f`
  n5=`sfget parform=n n5 < $f`
  n1=$((en1 - 88))
  n2=$((en2 - 88))
  n3=$((en3 - 88))
  f3=$((n3 / 2))
  < $f sfwindow f5=0 n5=${n5:-1} n4=1 f4=0 f3=$f3 n3=1 min1=0 min2=0 min3=0 n1=$n1 n2=$n2 | sfgrey color=j gainpanel=e title="$f" | sfpen
  #< $f sfwindow f5=0 n5=${n5:-1} n4=1 f4=0 f3=$f3 n3=1 min1=0 min2=0 min3=0 n1=$n1 n2=$n2 | sfattr
  #< $f sfwindow f5=0 n5=${n5:-1} n4=1 f4=0 f3=$f3 n3=1 min1=0 min2=0 min3=0 n1=$n1 n2=$n2 | sfgrey color=j gainpanel=e scalebar=y title="$f" > ${f/rsf/vpl}
}

show $1
