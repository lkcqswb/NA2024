import sys
import matplotlib.pyplot as plt

def read_points(file_name):
    points = []
    try:
        with open(file_name, 'r') as file:
            for line in file:
                point = list(map(float, line.strip().split()))
                points.append(point)
    except Exception as e:
        print(f"Error reading file: {e}")
    return points

def plot_points(points):
    if not points:
        print("No points to plot.")
        return

    # 将点分为 x 和 y 坐标
    x, y = zip(*points)

    # 绘制点
    plt.scatter(x, y, color='blue')
    plt.title('Scatter Plot of Points')
    plt.xlabel('X-axis')
    plt.ylabel('Y-axis')
    plt.grid()
    
    # 显示图形
    plt.show()

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python plot_points.py <filename>")
        sys.exit(1)

    file_name = sys.argv[1]
    points = read_points(file_name)
    plot_points(points)