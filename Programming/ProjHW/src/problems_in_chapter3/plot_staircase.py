import matplotlib.pyplot as plt
import numpy as np


def plot_subplots(data, order, min_val, max_val, file_name):
    # 创建一个 order x order 的子图网格
    fig, axs = plt.subplots(order, order, figsize=(10, 10))
    fig.subplots_adjust(hspace=0.4, wspace=0.4)

    # 遍历子图位置，控制下三角部分绘制
    for i in range(order):
        for j in range(order):
            if j <= i:  # 仅绘制下三角区域（包括对角线）
                axs[i, j].set_title(f't{i-1},{j}')
                x_vals = []
                y_vals = []
                for item in data:
                    i_1 = item[0] + 1
                    j_1 = item[1]
                    if i == i_1 and j == j_1:
                        x = item[2]
                        y = item[3]
                        x_vals.append(x)  # 获取 x 坐标
                        y_vals.append(y)  # 获取 y 坐标
                
                # 绘制这些点并将它们连接成线
                if x_vals and y_vals:
                    axs[i, j].plot(x_vals, y_vals, color='black', marker='o', linestyle='-', linewidth=1)
                
                # 设置x轴范围
                axs[i, j].set_xlim(min_val, max_val)  
                
                # 设置y轴范围：只使用当前子图的y值
                if y_vals:
                    axs[i, j].set_ylim(0, max(y_vals) * 1.1)  # 将y范围设置为最大y值的1.1倍
                axs[i, j].grid(True)  # 添加网格
            else:
                axs[i, j].axis('off')  # 关闭坐标轴，隐藏子图

    plt.savefig(file_name[:-4] + ".png")


def read_data(file_name):
    data = []
    order = 1
    min_val = 1e6
    max_val = -1e6
    import re
    pattern = re.compile(r"t_{\s*(-?\d+),\s*(-?\d+)\s*}(-?\d*\.?\d+(?:e-?\d+)?)\s*(-?\d*\.?\d+(?:e-?\d+)?)")
    with open(file_name, "r", encoding="utf-8") as file:
        for line in file:
            match = pattern.match(line.strip())  # 使用strip去除空白字符
            if match:
                x = int(match.group(1))
                if(order < x + 2): order = x + 2
                y = int(match.group(2))
                a = float(match.group(3))
                if(a > max_val): max_val = a
                if(a < min_val): min_val = a
                b = float(match.group(4))
                data.append((x, y, a, b))  
            else:
                print(f"未匹配的行: {line.strip()}")
    return data, order, min_val, max_val


if __name__ == "__main__":
    import sys
    if len(sys.argv) != 2:
        print("Usage: python plot_points.py <filename>")
        sys.exit(1)

    file_name = sys.argv[1]
    data, order, min_val, max_val = read_data(file_name)
    plot_subplots(data, order, min_val, max_val, file_name)
