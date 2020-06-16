for x in {ubuntu,fedora,debian}-{gcc,clang}-hepmc{2,3}-py{2,3}; do
    for t in "" "-latex"; do
        echo $x$t
        docker push hepstore/hepbase-$x$t
        echo
    done
done
