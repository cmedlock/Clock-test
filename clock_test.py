def parse_file(fname):
    # input: data file
    # output: x_copy,y_copy,x_command,y_command (np.arrays)
    x_copy,y_copy,x_command,y_command = [1,2,3],[1,2,3],[1,2,3],[1,2,3]
    return x_copy,y_copy,x_command,y_command
    
def plot_xyt_other(x,y,t,other,n,othername,fname):
    # input: x[t],y[t],t
    #        other[n],n (for e.g. |X[k]|,k)
    # output: xyt_othername_fname.png with plots of x[t] or y[t]
    #         directly above other[n]
    pass

def interpolate(a):
    # input: a (list) with missing points marked by the value -5
    # output: b (list) with the missing points filled in
    pass
    