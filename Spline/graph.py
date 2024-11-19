import numpy as np
import matplotlib.pyplot as plt

natural = np.loadtxt("natural_spline_error.txt")
true = np.loadtxt("true_spline_error.txt")

slope_natural, _ = np.polyfit(natural[:, 0], natural[:, 1], 1)
slope_true, _ = np.polyfit(true[:, 0], true[:, 1], 1)

print(f"Коэффициент наклона (Естественный сплайн): {slope_natural}")
print(f"Коэффициент наклона (Истинные границы): {slope_true}")

plt.plot(natural[:, 0], natural[:, 1], label="Естественный сплайн")
plt.plot(true[:, 0], true[:, 1], label="Истинные границы")
plt.legend()
plt.show()
