import matplotlib.pyplot as plt

''''# Pie chart
labels = ['LTR', 'DNA', 'SINE', 'LINE']
sizes = [8.84, 3.55,12.7,20.92]
# only "explode" the 2nd slice (i.e. 'Hogs')
explode = (0.05, 0.05, 0.05, 0.05)
fig1, ax1 = plt.subplots()
ax1.pie(sizes, explode=explode, autopct='%1.1f%%',
        shadow=True, startangle=90, labeldistance=1.15)
# Equal aspect ratio ensures that pie is drawn as a circle
ax1.axis('equal')
plt.tight_layout()
plt.savefig('my_svg_pie.svg', transparent=True)

plt.show()'''

plt.bar(['L','A','h'],[20,20,20])
plt.show()