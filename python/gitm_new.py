def read(filename, varlist=None, newfile=True):
    '''
    Read binary file from GITM output.
    '''

    from re import sub
    from struct import unpack
    import sys
    import datetime as dt
    import numpy as np

    out = {}
    print('Reading file...')
    if varlist is None:
        varlist = list()
    out['file'] = filename
    # Read data and header info
    f=open(filename, 'rb')

    # Using the first FORTRAN header, determine endian.
    # Default is little.
    out['endian']='little'
    endChar='>'
    rawRecLen=f.read(4)
    if len(rawRecLen) < 4:
        print("GitmBin ERROR: empty file [", filename, "]")
        sys.exit(1)
    recLen=(unpack(endChar+'l',rawRecLen))[0]
    if (recLen>10000)or(recLen<0):
        # Ridiculous record length implies wrong endian.
        out['endian']='big'
        endChar='<'
        recLen=(unpack(endChar+'l',rawRecLen))[0]

    # Read version; read fortran footer+header.
    out['version']=unpack(endChar+'d',f.read(recLen))[0]
    (oldLen, recLen)=unpack(endChar+'2l',f.read(8))

    # Read grid size information.
    (out['nLon'],out['nLat'],out['nAlt'])=\
        unpack(endChar+'lll',f.read(recLen))
    (oldLen, recLen)=unpack(endChar+'2l',f.read(8))

    # Read number of variables.
    out['nVars']=unpack(endChar+'l',f.read(recLen))[0]
    (oldLen, recLen)=unpack(endChar+'2l',f.read(8))

    # Collect variable names.
    var=[]
    for i in range(out['nVars']):
        var.append(unpack(endChar+'%is'%(recLen),f.read(recLen))[0])
        (oldLen, recLen)=unpack(endChar+'2l',f.read(8))

    # Extract time.
    (yy,mm,dd,hh,mn,ss,ms)=unpack(endChar+'lllllll',f.read(recLen))
    out['time']=dt.datetime(yy,mm,dd,hh,mn,ss,ms)
    (oldLen)=unpack(endChar+'l',f.read(4))


    # Read the rest of the data.
    nTotal=out['nLon']*out['nLat']*out['nAlt']
    var = [k.decode('utf-8') for k in var]
    for val in var:
        # Trim variable names.
        v=sub('\[|\]', '', val).strip()
        s=unpack(endChar+'l',f.read(4))[0]
        # Test to see if this variable is desired
        gvar=True
        if len(varlist) > 0:
            try:
                varlist.index(v)
            except ValueError:
                if((v.find('Altitude') < 0 and v.find('Longitude') < 0
                    and v.find('Latitude') < 0) or not newfile):
                    gvar=False
        # Unpack the data and save, if desired
        temp=unpack(endChar+'%id'%(nTotal),f.read(s))
        if gvar:
            out[v]=np.array(temp)
            # Reshape arrays, note that ordering in file is Fortran-like.
            out[v]=out[v].reshape(
                (out['nLon'],out['nLat'],out['nAlt']),
                order='fortran')

        f.read(4)
    out['dLat'] = out['Latitude']*180.0/np.pi
    out['dLon'] = out['Longitude']*180.0/np.pi
    ut = out['time'].hour+out['time'].minute/60+out['time'].second/3600
    out['LT'] = (ut + out['dLon']/15)%24
    return out
