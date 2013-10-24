from pylab import *
import os

datadir = 'data'
if __name__ == '__main__':
    data = {}
    for filename in os.listdir(datadir):
        path = os.path.join(datadir, filename)
        if os.path.isfile(path):
            if filename.startswith('GROUCHO'):
                celltype = filename.rpartition('.')[-1]
                data[celltype] = loadtxt(path)
                print 'Loaded', celltype
    field1mm_col = {'suppyrRS': 11,
                    'suppyrFRB': 10,
                    'tuftIB': 10,
                    'tuftRS': 10,
                    'nontuftRS': 10,
                }
    field2mm_col = {'suppyrRS': 12,
                'suppyrFRB': 11,
                'tuftIB': 11,
                'tuftRS': 11,
                'nontuftRS': 11,                
            }
    subplot(411)
    field_1mm = sum([data[celltype][:, field1mm_col[celltype]] for celltype in field1mm_col], axis=0)
    plot(data['suppyrFRB'][:,0], field_1mm, 'k-')
    title('Field at 1 mm')
    subplot(412)
    field_2mm = sum([data[celltype][:, field2mm_col[celltype]] for celltype in field2mm_col], axis=0)
    plot(data['suppyrFRB'][:,0], field_2mm, 'k-')
    title('Field at 2 mm')
    subplot(413)
    plot(data['suppyrRS'][:,0], data['suppyrRS'][:,1], 'k-', label='Sup. Pyr. RS')
    plot(data['suppyrFRB'][:,0], data['suppyrFRB'][:,1], 'r-', label='Sup. Pyr. FRB')
    subplot(414)
    plot(data['spinstell'][:,0], data['spinstell'][:,1], 'k-', label='Spiny Stellate')
    plot(data['tuftIB'][:,0], data['tuftIB'][:,1], 'r-', label='tufted IB')
    legend(loc='best')
    show()
