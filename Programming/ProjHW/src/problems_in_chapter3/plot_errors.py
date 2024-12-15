# -*- coding: utf-8 -*-  # 在文件头部添加编码声明

import matplotlib.pyplot as plt
import csv

# 读取 CSV 文件
def read_errors(filename="errors.csv"):
    Ns = []
    errors_pp = []
    errors_bs = []
    with open(filename, mode='r', encoding='utf-8') as file:  # 确保读取时使用 UTF-8 编码
        csv_reader = csv.reader(file)
        next(csv_reader)  # 跳过表头
        for row in csv_reader:
            Ns.append(int(row[0]))
            errors_pp.append(float(row[1]))
            errors_bs.append(float(row[2]))
    return Ns, errors_pp, errors_bs

# 绘制最大误差图
def plot_errors():
    Ns, errors_pp, errors_bs = read_errors()

    # 绘制误差
    plt.plot(Ns, errors_pp, label='PPForm Error', marker='o', linestyle='-', color='b')
    plt.plot(Ns, errors_bs, label='BSpline Error', marker='o', linestyle='--', color='r')

    # 添加图例和标签
    plt.xlabel('Number of Points (N)')
    plt.ylabel('Max Error')
    plt.title('Max Error vs. Number of Points')
    plt.legend()

    # 显示图形
    plt.grid(True)
    plt.savefig("max_error_plot.png")  # 将图形保存为文件

# 调用函数绘制图形
plot_errors()
