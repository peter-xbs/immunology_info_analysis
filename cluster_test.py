import matplotlib.pyplot as plt
import seaborn as sns
# sns.set(color_codes=True)
# iris = sns.load_dataset("iris")
# species = iris.pop("species")
# #设置图片大小
# g= sns.clustermap(iris, fmt="d",cmap='YlGnBu',figsize=(6,9))
# ax = g.ax_heatmap
# label_y = ax.get_yticklabels(minor=True)
# plt.setp(label_y, rotation=360, horizontalalignment='left')
# #设置图片名称，分辨率，并保存
# plt.savefig('cluster.pdf',dpi = 300)
# plt.show()

# import math, datetime
# import rpy2.robjects.lib.ggplot2 as ggplot2
# import rpy2.robjects as ro
# from rpy2.robjects.packages import importr
# base = importr('base')
# from rpy2.robjects import r
#
#
# mtcars = data(datasets).fetch('mtcars')['mtcars']
# pp = ggplot2.ggplot(mtcars) + \
#      ggplot2.aes_string(x='wt', y='mpg', col='factor(cyl)') + \
#      ggplot2.geom_point() + \
#      ggplot2.geom_smooth(ggplot2.aes_string(group = 'cyl'),
#                          method = 'lm')
# pp.plot()

