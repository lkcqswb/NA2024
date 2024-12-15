import matplotlib.pyplot as plt
import sys
# 读取 Pc.txt 文件
def read_data(filename):
    x_vals, y_vals = [], []
    labels = []  # 用于存储每段数据的标签
    with open(filename, 'r') as file:
        lines = file.readlines()
        temp_x, temp_y = [], []
        for line in lines:
            if line.startswith("#END#"):
                if temp_x:
                    x_vals.append(temp_x)
                    y_vals.append(temp_y)
                    labels.append(line[len('#END#'):].strip())
                temp_x, temp_y = [], []  # 清空临时列表，准备下一段数据
            elif not line.startswith("#"):
                x, y = map(float, line.split())
                temp_x.append(x)
                temp_y.append(y)
        if temp_x:  # 处理最后一段数据
            x_vals.append(temp_x)
            y_vals.append(temp_y)
    return x_vals, y_vals, labels

# 绘制图形
def plot_data(file_name):
    x_vals, y_vals, labels = read_data(file_name)
    for i in range(len(x_vals)):
        plt.plot(x_vals[i], y_vals[i], label=labels[i], marker='o')
    
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.title('BSpline and Function Comparison')
    plt.legend()
    plt.grid(True)
    plt.savefig(file_name[:-4]+".png")


if __name__ == "__main__":
    
    if len(sys.argv) != 2:
        print("Usage: python plot_points.py <filename>")
        sys.exit(1)

    file_name = sys.argv[1]
    plot_data(file_name)
