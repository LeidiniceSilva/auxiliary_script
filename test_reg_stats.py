import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as s
import matplotlib.lines as mlines
import matplotlib.transforms as mtransforms

from scipy.stats import gamma

data_a = [[1,2,5], [5,7,2,2,5], [7,2,5]]
data_b = [[6,4,2], [1,2,5,3,2], [2,3,5,1]]
data_c = [[6,4,2], [1,2,5,3,2], [2,3,5,1]]
data_d = [[6,4,2], [1,2,5,3,2], [2,3,5,1]]

ticks = ['Historical', 'RCP2.6', 'RCP8.5']

def set_box_color(bp, color):
    plt.setp(bp['boxes'], color=color)
    plt.setp(bp['whiskers'], color=color)
    plt.setp(bp['caps'], color=color)
    plt.setp(bp['medians'], color=color)

ax1 = plt.subplot(3, 1, 1)
plt.axhline(2, linewidth=1., linestyle='dashed', color='black')
plt.axhline(4, linewidth=1., linestyle='dashed', color='black')
plt.axhline(6, linewidth=1., linestyle='dashed', color='black')
plt.axhline(8, linewidth=1., linestyle='dashed', color='black')
bpa = plt.boxplot(data_a, positions=np.array(range(len(data_a)))*2.0-0.55, sym='', widths=0.25)
bpb = plt.boxplot(data_b, positions=np.array(range(len(data_b)))*2.0-0.15, sym='', widths=0.25)
bpc = plt.boxplot(data_c, positions=np.array(range(len(data_c)))*2.0+0.15, sym='', widths=0.25)
bpd = plt.boxplot(data_d, positions=np.array(range(len(data_d)))*2.0+0.55, sym='', widths=0.25)
set_box_color(bpa, 'blue') 
set_box_color(bpb, 'red')
set_box_color(bpc, 'green')
set_box_color(bpd, 'gray')

plt.xlim(-2, len(ticks)*2)
plt.xticks(range(0, len(ticks) * 2, 2), ticks, fontweight='bold')
plt.ylim(0, 10)
plt.yticks(np.arange(0, 12, 2))
ax1.set_xticklabels([])

ax2 = plt.subplot(3, 1, 2)
plt.axhline(2, linewidth=1., linestyle='dashed', color='black')
plt.axhline(4, linewidth=1., linestyle='dashed', color='black')
plt.axhline(6, linewidth=1., linestyle='dashed', color='black')
plt.axhline(8, linewidth=1., linestyle='dashed', color='black')
bpa = plt.boxplot(data_a, positions=np.array(range(len(data_a)))*2.0-0.55, sym='', widths=0.25)
bpb = plt.boxplot(data_b, positions=np.array(range(len(data_b)))*2.0-0.15, sym='', widths=0.25)
bpc = plt.boxplot(data_c, positions=np.array(range(len(data_c)))*2.0+0.15, sym='', widths=0.25)
bpd = plt.boxplot(data_d, positions=np.array(range(len(data_d)))*2.0+0.55, sym='', widths=0.25)
set_box_color(bpa, 'blue') 
set_box_color(bpb, 'red')
set_box_color(bpc, 'green')
set_box_color(bpd, 'gray')

plt.ylabel(u'Precipitação (mm d⁻¹)', fontweight='bold')
plt.xlim(-2, len(ticks)*2)
plt.xticks(range(0, len(ticks) * 2, 2), ticks, fontweight='bold')
plt.ylim(0, 10)
plt.yticks(np.arange(0, 12, 2))
ax2.set_xticklabels([])

ax3 = plt.subplot(3, 1, 3)
plt.axhline(2, linewidth=1., linestyle='dashed', color='black')
plt.axhline(4, linewidth=1., linestyle='dashed', color='black')
plt.axhline(6, linewidth=1., linestyle='dashed', color='black')
plt.axhline(8, linewidth=1., linestyle='dashed', color='black')
bpa = plt.boxplot(data_a, positions=np.array(range(len(data_a)))*2.0-0.55, sym='', widths=0.25)
bpb = plt.boxplot(data_b, positions=np.array(range(len(data_b)))*2.0-0.15, sym='', widths=0.25)
bpc = plt.boxplot(data_c, positions=np.array(range(len(data_c)))*2.0+0.15, sym='', widths=0.25)
bpd = plt.boxplot(data_d, positions=np.array(range(len(data_d)))*2.0+0.55, sym='', widths=0.25)
set_box_color(bpa, 'blue') 
set_box_color(bpb, 'red')
set_box_color(bpc, 'green')
set_box_color(bpd, 'gray')

plt.xlim(-2, len(ticks)*2)
plt.xticks(range(0, len(ticks) * 2, 2), ticks, fontweight='bold')
plt.ylim(0, 10)
plt.yticks(np.arange(0, 12, 2))
plt.show()
exit()

fig=plt.figure()

ax1 = plt.subplot(3, 5, 1)
plt.scatter(3, 10, s=80, c='blue', marker='o')
plt.scatter(3, 14, s=80, c='blue', marker='D')
plt.scatter(-2, -10, s=80, c='red', marker='o')
plt.scatter(1, -5, s=80, c='red', marker='D')
ax1.set_xticklabels([])
plt.axvline(0, linewidth=1., linestyle='dashed', color='black')
plt.axhline(0, linewidth=1., linestyle='dashed', color='black')
title = ax1.text(-2.0,17.3, '   Annual  ', fontweight='bold', zorder=6, color='k',
                bbox={'facecolor':'silver', 'alpha':0.5, 'pad':4})
             
ax2 = plt.subplot(3, 5, 2)
plt.scatter(3, 10, s=80, c='blue', marker='o')
plt.scatter(3, 14, s=80, c='blue', marker='D')
plt.scatter(-2, -10, s=80, c='red', marker='o')
plt.scatter(1, -5, s=80, c='red', marker='D')
ax2.set_xticklabels([])
ax2.set_yticklabels([])
plt.axvline(0, linewidth=1., linestyle='dashed', color='black')
plt.axhline(0, linewidth=1., linestyle='dashed', color='black')
title = ax2.text(-2.0,17.3, '     DJF     ', fontweight='bold', zorder=6, color='k',
                bbox={'facecolor':'silver', 'alpha':0.5, 'pad':4})
                
ax3 = plt.subplot(3, 5, 3)
plt.scatter(3, 10, s=80, c='blue', marker='o')
plt.scatter(3, 14, s=80, c='blue', marker='D')
plt.scatter(-2, -10, s=80, c='red', marker='o')
plt.scatter(1, -5, s=80, c='red', marker='D')
ax3.set_xticklabels([])
ax3.set_yticklabels([])
plt.axvline(0, linewidth=1., linestyle='dashed', color='black')
plt.axhline(0, linewidth=1., linestyle='dashed', color='black')
title = ax3.text(-2.0,17.3, '    MAM    ', fontweight='bold', zorder=6, color='k',
                bbox={'facecolor':'silver', 'alpha':0.5, 'pad':4})

ax4 = plt.subplot(3, 5, 4)
plt.scatter(3, 10, s=80, c='blue', marker='o')
plt.scatter(3, 14, s=80, c='blue', marker='D')
plt.scatter(-2, -10, s=80, c='red', marker='o')
plt.scatter(1, -5, s=80, c='red', marker='D')
ax4.set_xticklabels([])
ax4.set_yticklabels([])
plt.axvline(0, linewidth=1., linestyle='dashed', color='black')
plt.axhline(0, linewidth=1., linestyle='dashed', color='black')
title = ax4.text(-2.0,17.3, '      JJA     ', fontweight='bold', zorder=6, color='k',
                bbox={'facecolor':'silver', 'alpha':0.5, 'pad':4})
   
ax5 = plt.subplot(3, 5, 5)
plt.scatter(3, 10, s=80, c='blue', marker='o', label='Had (RCP2.6)')
plt.scatter(3, 14, s=80, c='blue', marker='D', label='Reg_Had (RCP2.6)')
plt.scatter(-2, -10, s=80, c='red', marker='o', label='Had (RCP8.55)')
plt.scatter(1, -5, s=80, c='red', marker='D', label='Reg_Had (RCP8.55)')
ax5.set_xticklabels([])
ax5.set_yticklabels([])
plt.axvline(0, linewidth=1., linestyle='dashed', color='black')
plt.axhline(0, linewidth=1., linestyle='dashed', color='black')
title = ax5.text(-2.0,17.3, '     SON    ', fontweight='bold', zorder=6, color='k',
                bbox={'facecolor':'silver', 'alpha':0.5, 'pad':4})
title = ax5.text(3.75,11., '      AMZ       ', fontweight='bold', rotation=270, zorder=6, color='k',
                bbox={'facecolor':'silver', 'alpha':0.5, 'pad':4})
      
ax6 = plt.subplot(3, 5, 6)
plt.scatter(3, 10, s=80, c='blue', marker='o')
plt.scatter(3, 14, s=80, c='blue', marker='D')
plt.scatter(-2, -10, s=80, c='red', marker='o')
plt.scatter(1, -5, s=80, c='red', marker='D')
ax6.set_xticklabels([])
plt.axvline(0, linewidth=1., linestyle='dashed', color='black')
plt.axhline(0, linewidth=1., linestyle='dashed', color='black')
plt.ylabel(u'Precipitation change (%)', fontweight='bold')

ax7 = plt.subplot(3, 5, 7)
plt.scatter(3, 10, s=80, c='blue', marker='o')
plt.scatter(3, 14, s=80, c='blue', marker='D')
plt.scatter(-2, -10, s=80, c='red', marker='o')
plt.scatter(1, -5, s=80, c='red', marker='D')
ax7.set_xticklabels([])
ax7.set_yticklabels([])
plt.axvline(0, linewidth=1., linestyle='dashed', color='black')
plt.axhline(0, linewidth=1., linestyle='dashed', color='black')

ax8 = plt.subplot(3, 5, 8)
plt.scatter(3, 10, s=80, c='blue', marker='o')
plt.scatter(3, 14, s=80, c='blue', marker='D')
plt.scatter(-2, -10, s=80, c='red', marker='o')
plt.scatter(1, -5, s=80, c='red', marker='D')
ax8.set_xticklabels([])
ax8.set_yticklabels([])
plt.axvline(0, linewidth=1., linestyle='dashed', color='black')
plt.axhline(0, linewidth=1., linestyle='dashed', color='black')

ax9 = plt.subplot(3, 5, 9)
plt.scatter(3, 10, s=80, c='blue', marker='o')
plt.scatter(3, 14, s=80, c='blue', marker='D')
plt.scatter(-1, -10, s=80, c='red', marker='o')
plt.scatter(1, -5, s=80, c='red', marker='D')
ax9.set_xticklabels([])
ax9.set_yticklabels([])
plt.axvline(0, linewidth=1., linestyle='dashed', color='black')
plt.axhline(0, linewidth=1., linestyle='dashed', color='black')

ax10 = plt.subplot(3, 5, 10)
plt.scatter(3, 10, s=80, c='blue', marker='o')
plt.scatter(3, 14, s=80, c='blue', marker='D')
plt.scatter(-2, -10, s=80, c='red', marker='o')
plt.scatter(1, -5, s=80, c='red', marker='D')
ax10.set_xticklabels([])
ax10.set_yticklabels([])
plt.axvline(0, linewidth=1., linestyle='dashed', color='black')
plt.axhline(0, linewidth=1., linestyle='dashed', color='black')
title = ax10.text(3.75,11., '      NEB       ', fontweight='bold', rotation=270, zorder=6, color='k',
                bbox={'facecolor':'silver', 'alpha':0.5, 'pad':4})
      
ax11 = plt.subplot(3, 5, 11)
plt.scatter(3, 10, s=80, c='blue', marker='o')
plt.scatter(3, 14, s=80, c='blue', marker='D')
plt.scatter(-2, -10, s=80, c='red', marker='o')
plt.scatter(1, -5, s=80, c='red', marker='D')
plt.axvline(0, linewidth=1., linestyle='dashed', color='black')
plt.axhline(0, linewidth=1., linestyle='dashed', color='black')

ax12 = plt.subplot(3, 5, 12)
plt.scatter(3, 10, s=80, c='blue', marker='o')
plt.scatter(3, 14, s=80, c='blue', marker='D')
plt.scatter(-2, -10, s=80, c='red', marker='o')
plt.scatter(1, -5, s=80, c='red', marker='D')
ax12.set_yticklabels([])
plt.axvline(0, linewidth=1., linestyle='dashed', color='black')
plt.axhline(0, linewidth=1., linestyle='dashed', color='black')

ax13 = plt.subplot(3, 5, 13)
plt.scatter(3, 10, s=80, c='blue', marker='o')
plt.scatter(3, 14, s=80, c='blue', marker='D')
plt.scatter(-2, -10, s=80, c='red', marker='o')
plt.scatter(1, -5, s=80, c='red', marker='D')
ax13.set_yticklabels([])
plt.axvline(0, linewidth=1., linestyle='dashed', color='black')
plt.axhline(0, linewidth=1., linestyle='dashed', color='black')
plt.xlabel(u'Temperature change (°C)', fontweight='bold')

ax14 = plt.subplot(3, 5, 14)
plt.scatter(3, 10, s=80, c='blue', marker='o')
plt.scatter(3, 14, s=80, c='blue', marker='D')
plt.scatter(-2, -10, s=80, c='red', marker='o')
plt.scatter(1, -5, s=80, c='red', marker='D')
ax14.set_yticklabels([])
plt.axvline(0, linewidth=1., linestyle='dashed', color='black')
plt.axhline(0, linewidth=1., linestyle='dashed', color='black')

ax15 = plt.subplot(3, 5, 15)
plt.scatter(3, 10, s=80, c='blue', marker='o')
plt.scatter(3, 14, s=80, c='blue', marker='D')
plt.scatter(-2, -10, s=80, c='red', marker='o')
plt.scatter(1, -5, s=80, c='red', marker='D')
ax15.set_yticklabels([])
plt.axvline(0, linewidth=1., linestyle='dashed', color='black')
plt.axhline(0, linewidth=1., linestyle='dashed', color='black')
title = ax15.text(3.75,11., ' MATOPIBA  ', fontweight='bold', rotation=270, zorder=6, color='k',
                bbox={'facecolor':'silver', 'alpha':0.6, 'pad':4})
                
ax5.legend(bbox_to_anchor=(1.4, 1), loc=2, borderaxespad=0.5)
plt.show()
exit()

# dataset-1 
x1 = [3, 2, 3, 3, 4, 2, 4, 2, 1, 2] 
y1 = [3, 3, 3, 3, 4, 2, 3, 1, 1, 2] 
  
# dataset2 
x2 = [7, 6, 6, 7, 7, 7, 6, 5, 5, 5] 
y2 = [7, 7, 6, 7, 7, 7, 6, 5, 5, 5] 
 
# dataset2 
x3 = [8, 9, 9, 9, 10, 9, 9, 8, 9, 8] 
y3 = [8, 8, 9, 10, 8, 8, 9, 8, 10, 9] 

fig, ax = plt.subplots()
line = mlines.Line2D([0, 1], [0, 1], lw=1., color='black')
transform = ax.transAxes
line.set_transform(transform)
ax.add_line(line)  

plt.scatter(x1, y1, c ="red",  
            linewidths = 2.,  
            marker ="o",  
            edgecolor ="gray",  
            s = 50.) 
  
plt.scatter(x2, y2, c ="blue", 
            linewidths = 2., 
            marker ="o",  
            edgecolor ="gray",  
            s = 50.) 

plt.scatter(x3, y3, c ="green", 
            linewidths = 2., 
            marker ="o",  
            edgecolor ="gray",  
            s = 50.) 
            
plt.xlabel("X-axis") 
plt.ylabel("Y-axis") 
plt.grid(True, lw=1., ls ='--', color='gray')
plt.show()
exit()

fig, ax = plt.subplots(1, 1)

# Calculate a few first moments:
a = 1.99323054838
mean, var, skew, kurt = gamma.stats(a, moments='mvsk')

# Display the probability density function (pdf):
x = np.linspace(gamma.ppf(0.01, a), gamma.ppf(0.99, a), 100)
ax.plot(x, gamma.pdf(x, a), 'r-', lw=5, alpha=0.6, label='gamma pdf')

# Alternatively, freeze the distribution and display the frozen pdf:
rv = gamma(a)
ax.plot(x, rv.pdf(x), 'k-', lw=2, label='frozen pdf')
r = gamma.rvs(a, size=1000)

# And compare the histogram:
ax.hist(r, normed=True, histtype='stepfilled', color='gray', alpha=0.8)
ax.legend(loc='best', frameon=False)

plt.show()
exit()

wei = s.weibull_min(2, 0, 2) # shape, loc, scale - creates weibull object
sample = wei.rvs(1000)
shape, loc, scale = s.weibull_min.fit(sample, floc=0) 
x = np.linspace(np.min(sample), np.max(sample))

plt.hist(sample, normed=True, fc="none", ec="grey", label="frequency")
plt.plot(x, wei.cdf(x), label="cdf")
plt.plot(x, wei.pdf(x), label="pdf")
plt.plot(x[1:], np.diff(wei.cdf(x)), label="derivative")
plt.legend(loc=1)
plt.show()
