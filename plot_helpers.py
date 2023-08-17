import config


def hover_points(points, threshold=config.THRESHOLD):
    ### Return points to interact with via hover if there are fewer than threshold
    if len(points) > threshold:
        return points.iloc[:0]
    return points


def filter_points(points, x_range, y_range):
    if x_range and y_range:
        return points[x_range, y_range]
    return points