1.include中定义了位错mobility law，采用BCC W热激活滑移模型
2.在sizeLLx 文件夹下可以单独编译、运行,matlab运行viz.m可视化
    make: 编译当前文件夹
    ./dd_sim: 运行当前文件夹
    make clean: 清除编译文件
    make empty: 删除当前文件夹下的运行结果
3.在bending_sizes 文件夹下对所有子文件夹进行编译，matlab运行viz4allsize.m可视化
    make build: 编译所有子文件夹
    make run: 运行所有子文件夹
    make clean: 清除编译文件
    make empty: 删除所有子文件夹下的运行结果
4.input.txt: 输入文件