import csv
import os


def write_to_csvfile(filename, data):

    file = f'{filename}.csv'
    if not os.path.isfile(file):
        write_file = open(file, 'w+')
        writer = csv.writer(write_file, delimiter=',')
        writer.writerow(['Cores', 'Resolution', 'Epoch', 'Time'])
        write_file.close()

    with open(file, 'a') as write_file:
        writer = csv.writer(write_file)
        writer.writerow(data)
