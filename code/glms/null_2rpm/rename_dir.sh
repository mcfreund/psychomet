files=($(find -name *_null_2rpm_))
for file_i in ${!files[@]}; do
    mv ${files[$file_i]} ${files[$file_i]::-1}
done
