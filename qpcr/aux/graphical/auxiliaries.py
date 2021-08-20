import plotly.express as px
import plotly as py
import matplotlib.pyplot as plt


#================== Layout Auxiliaries ==================
def adjust_layout(graph_number):
    if graph_number%2 == 0:
            if graph_number%4 == 0:
                rows = int(graph_number/4)
                cols = int(graph_number/rows)
            else:
                rows = int(graph_number/2)
                cols = int(graph_number/rows)
    elif graph_number%3 == 0:
        rows = int(graph_number/3)
        cols = int(graph_number/rows)
    elif graph_number%5 == 0:
        rows = graph_number/5
        cols = int(graph_number/rows)
    else:
        rows = graph_number
        cols = 1
    return int(rows), int(cols)