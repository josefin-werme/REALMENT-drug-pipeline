# command used to compress/archive files
# for I in $(cat ../tissues.txt); do echo $I; tar -I "zstd -T0" -cf $I.tar.zst $I; done

# to extract files, use
# for I in $(cat ../tissues.txt); do echo $I; tar -I "zstd -T0" -xf $I.tar.zst $I; done
