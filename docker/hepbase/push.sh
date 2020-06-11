for x in {ubuntu,fedora,debian}-{gcc,clang}-hepmc{2,3}-py{2,3}; do
    echo $x
    docker push hepstore/hepbase-$x
    echo
done
