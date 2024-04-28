import time


def calculate_time(func):

    def wrapper(*args, **kwargs):
        start_time = time.time()
        func(*args, **kwargs)
        end_time = time.time()
        print(f"Time taken is: {end_time - start_time}")

    return wrapper
