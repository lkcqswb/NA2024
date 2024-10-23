import matplotlib.pyplot as plt
import sys

def plot_data(x_values, y_values, label):
    plt.plot(x_values, y_values, label=label)

def plot_point(x_values, y_values, label):
    plt.scatter(x_values, y_values, label=label)

if len(sys.argv) != 2:
    print("Usage: python plot.py <filename>")
    sys.exit(1)

filename = sys.argv[1]

with open(filename, 'r') as file:
    lines = file.readlines()

curve_label = ""
x_values, y_values = [], []

for line in lines:
    line = line.strip()
    if line.startswith('#END#dots'):
        curve_label = line[len('#END#'):].strip()
        if x_values and y_values:
            plot_point(x_values, y_values, " characteristic points.")
        x_values, y_values = [], []
    elif line.startswith('#END#'):
        curve_label = line[len('#END#'):].strip()
        if x_values and y_values:
            plot_data(x_values, y_values, curve_label)
        x_values, y_values = [], []
    else:
        if line: 
            parts = line.split()
            if len(parts) == 2:
                x, y = map(float, parts)
                x_values.append(x)
                y_values.append(y)



plt.xlabel('x')
plt.ylabel('y')
plt.title(f"{filename[:-4]}_"+curve_label+"_Plots")
plt.legend()
plt.grid(True)
plt.savefig(f"{filename[:-4]}"+"_"+curve_label+"_"+"plot.png")
print("The image has been saved to the current directory.")