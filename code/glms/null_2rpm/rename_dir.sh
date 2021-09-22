files=($(find -name *_null_2rpm_))
for file_i in ${!files[@]}; do
    mv ${files[$file_i]} ${files[$file_i]::-1}
done


files=($(find out/glms -maxdepth 4 -type d -name *null_2rpm_wave2))
echo ${files[@]}
for file_i in ${!files[@]}; do
    name_file=${files[$file_i]}
    name_file_edit=${files[$file_i]::-1}
    name_file_new=${name_file_edit}1
    #echo $name_file_new
    mv $name_file $name_file_new
done

unset name_file name_file_edit name_file_new files

files=($(find out/glms -maxdepth 4 -type d -name *null_2rpm_wave3))
echo ${files[@]}
for file_i in ${!files[@]}; do
    name_file=${files[$file_i]}
    name_file_edit=${files[$file_i]::-1}
    name_file_new=${name_file_edit}2
    #echo $name_file_new
    mv $name_file $name_file_new
done

