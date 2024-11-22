import matplotlib.pyplot as plt

# Sample data
labels = ['Category A', 'Category B', 'Category C', 'Category D']
sizes = [25, 30, 20, 25]

# Colors for each category
colors = ['gold', 'lightcoral', 'lightskyblue', 'lightgreen']

# Exploding a slice (optional)
explode = (0.1, 0.1, 0.1, 0.1)

# Plotting the pie chart
plt.pie(sizes, explode=explode,  autopct='%1.1f%%', shadow=True, startangle=140)

# Equal aspect ratio ensures that pie is drawn as a circle
plt.axis('equal')

# Adding a title
plt.title('Example Pie Chart')

# Display the plot
plt.show()