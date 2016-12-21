import numpy as np
import matplotlib.pyplot as plt

col_labels = [  'Nome', 'Coluna 2', 'Coluna 3']
table_vals = [ [  'C1',         11,         12], 
               [  'C2',         21,         22], 
               [  'C3',         31,         32], 
               [  'C4',         41,         42], 
               [  'C5',         51,         52] ]

fig=plt.figure(figsize=(7,3.5))
axs = fig.add_subplot(111)
axs.axis('tight')
axs.axis('off')

the_table = axs.table(cellText=table_vals, colLoc='center',
                      colLabels=col_labels, 
                      colWidths=[.4,.2,.2], 
                      bbox=[-0.15,-0.1,1.27,1.2])  
                      
# bbox=[left, bottom, width, height]



plt.savefig("table.png")
# plt.show()
