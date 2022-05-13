# define EFS, read and write EFS
# Modified from codes of Daniel Trugman and Robin Matoza
# Wenyuan Fan 2021, 03, 01
#
# IV modifications:
# 
# switched sign of field tshead['tdif']
# changed filed tshead['deldist'] from deg to km to match documentation 


import numpy as np
import matplotlib.pyplot as plt
import struct
import obspy


#### EFS Class Definition ###
class EFS():
    '''
    Class definition for EFS-format data.
    Basic initialization syntax: edata = EFS("path_to_efs_file").
    '''

    ### --- Initialization --- ###
    #  efsfname: name of file to read
    #  prec_wf: precision for waveform arrays, 32 or 64
    #  prec_bp: precision for byteposition arrays
    def __init__(self, efsfname=None, prec_wf=np.float32, prec_bp=np.int32):

        # initialize fields
        self.fhead = {}
        self.ehead = {}
        self.waveforms = []

        # return here without file
        if efsfname is None:
            return

        # Open the EFS binary file
        f = open(efsfname, 'rb')

        # Assemble file header (20 bytes)
        self.fhead['bytetype'] = struct.unpack('i', f.read(4))[0]
        self.fhead['eheadtype'] = struct.unpack('i', f.read(4))[0]
        self.fhead['nbytes_ehead'] = struct.unpack('i', f.read(4))[0]
        self.fhead['tsheadtype'] = struct.unpack('i', f.read(4))[0]
        self.fhead['nbytes_tshead'] = struct.unpack('i', f.read(4))[0]

        # Assemble event header (184 bytes)
        self.ehead['efslabel'] = f.read(40).decode('UTF-8')
        self.ehead['datasource'] = f.read(40).decode('UTF-8')
        self.ehead['maxnumts'] = struct.unpack('i', f.read(4))[0]
        self.ehead['numts'] = struct.unpack('i', f.read(4))[0]
        self.ehead['cuspid'] = struct.unpack('i', f.read(4))[0]
        self.ehead['qtype'] = f.read(4).decode('UTF-8')
        self.ehead['qmag1type'] = f.read(4).decode('UTF-8')
        self.ehead['qmag2type'] = f.read(4).decode('UTF-8')
        self.ehead['qmag3type'] = f.read(4).decode('UTF-8')
        self.ehead['qmomenttype'] = f.read(4).decode('UTF-8')
        self.ehead['qlocqual'] = f.read(4).decode('UTF-8')
        self.ehead['qfocalqual'] = f.read(4).decode('UTF-8')
        self.ehead['qlat'] = struct.unpack('f', f.read(4))[0]
        self.ehead['qlon'] = struct.unpack('f', f.read(4))[0]
        self.ehead['qdep'] = struct.unpack('f', f.read(4))[0]
        self.ehead['qsc'] = struct.unpack('f', f.read(4))[0]
        self.ehead['qmag1'] = struct.unpack('f', f.read(4))[0]
        self.ehead['qmag2'] = struct.unpack('f', f.read(4))[0]
        self.ehead['qmag3'] = struct.unpack('f', f.read(4))[0]
        self.ehead['qmoment'] = struct.unpack('f', f.read(4))[0]
        self.ehead['qstrike'] = struct.unpack('f', f.read(4))[0]
        self.ehead['qdip'] = struct.unpack('f', f.read(4))[0]
        self.ehead['qrake'] = struct.unpack('f', f.read(4))[0]
        self.ehead['qyr'] = struct.unpack('i', f.read(4))[0]
        self.ehead['qmon'] = struct.unpack('i', f.read(4))[0]
        self.ehead['qdy'] = struct.unpack('i', f.read(4))[0]
        self.ehead['qhr'] = struct.unpack('i', f.read(4))[0]
        self.ehead['qmn'] = struct.unpack('i', f.read(4))[0]

        # 20 4-byte fields reserved for future uses - skip (80 bytes)
        for idum in range(0, 20):
            dummy = struct.unpack('i', f.read(4))[0]

        # Get byte positions for all time series (i64 or i32 arrays)
        bytepos = np.fromfile(f, dtype = prec_bp, count = self.ehead['numts'])
        self.ehead['bytepos'] = bytepos

        # Now loop over all the time series
        for ii in range(0, len(bytepos)):

            # Assemble tshead
            f.seek(bytepos[ii])
            # print(bytepos[ii])
            tshead = {}
            tshead['stname'] = f.read(8).decode('UTF-8')
            # print(tshead['stname'])
            tshead['loccode'] = f.read(8).decode('UTF-8')
            tshead['datasource'] = f.read(8).decode('UTF-8')
            tshead['sensor'] = f.read(8).decode('UTF-8')
            tshead['units'] = f.read(8).decode('UTF-8')
            tshead['chnm'] = f.read(4).decode('UTF-8')
            tshead['stype'] = f.read(4).decode('UTF-8')
            tshead['dva'] = f.read(4).decode('UTF-8')
            tshead['pick1q'] = f.read(4).decode('UTF-8')
            tshead['pick2q'] = f.read(4).decode('UTF-8')
            tshead['pick3q'] = f.read(4).decode('UTF-8')
            tshead['pick4q'] = f.read(4).decode('UTF-8')
            tshead['pick1name'] = f.read(4).decode('UTF-8')
            tshead['pick2name'] = f.read(4).decode('UTF-8')
            tshead['pick3name'] = f.read(4).decode('UTF-8')
            tshead['pick4name'] = f.read(4).decode('UTF-8')
            tshead['ppolarity'] = f.read(4).decode('UTF-8')
            tshead['problem'] = f.read(4).decode('UTF-8')
            # print(tshead)
            tshead['npts'] = struct.unpack('i', f.read(4))[0]
            tshead['syr'] = struct.unpack('i', f.read(4))[0]
            tshead['smon'] = struct.unpack('i', f.read(4))[0]
            tshead['sdy'] = struct.unpack('i', f.read(4))[0]
            tshead['shr'] = struct.unpack('i', f.read(4))[0]
            tshead['smn'] = struct.unpack('i', f.read(4))[0]
            tshead['compazi'] = struct.unpack('f', f.read(4))[0]
            tshead['compang'] = struct.unpack('f', f.read(4))[0]
            tshead['gain'] = struct.unpack('f', f.read(4))[0]
            tshead['f1'] = struct.unpack('f', f.read(4))[0]
            tshead['f2'] = struct.unpack('f', f.read(4))[0]
            tshead['dt'] = struct.unpack('f', f.read(4))[0]
            tshead['ssc'] = struct.unpack('f', f.read(4))[0]
            tshead['tdif'] = struct.unpack('f', f.read(4))[0]
            tshead['slat'] = struct.unpack('f', f.read(4))[0]
            tshead['slon'] = struct.unpack('f', f.read(4))[0]
            tshead['selev'] = struct.unpack('f', f.read(4))[0]
            tshead['deldist'] = struct.unpack('f', f.read(4))[0]
            tshead['sazi'] = struct.unpack('f', f.read(4))[0]
            tshead['qazi'] = struct.unpack('f', f.read(4))[0]
            tshead['pick1'] = struct.unpack('f', f.read(4))[0]
            tshead['pick2'] = struct.unpack('f', f.read(4))[0]
            tshead['pick3'] = struct.unpack('f', f.read(4))[0]
            tshead['pick4'] = struct.unpack('f', f.read(4))[0]

            # 20 4-byte fields reserved for future uses - skip
            for idum in range(0, 20):
                dummy = struct.unpack('i', f.read(4))[0]

            # Read the time-series itself
            data = np.fromfile(f, dtype = prec_wf, count = tshead['npts'])

            # Bundle tsheader and time-series for this waveform into efsdata, then append to list
            efsdata = tshead
            efsdata['data'] = data
            self.waveforms.append(efsdata)

    ###########################################################

    ### Function to convert from EFS to ObsPy Stream
    def to_obspy(self, keep_evdata=True, keep_stdata=True, keep_pkdata=True):

        '''
        Function to return obspy stream from EFS.
        Optional arguments preserve event, station, pick information in the
        stats dictionary for each trace.
        '''

        # Error checking
        try:
            from obspy.core import Stream, Trace, UTCDateTime
        except:
            print('EFS.to_obspy() is not available. ObsPy is not installed.')
            raise

        # Initialize stream
        st = Stream()

        # Optional: preserve event information in header
        if keep_evdata:
            evdata = {}
            evdata['evid'] = self.ehead['cuspid']
            for key in ['qlat', 'qlon', 'qdep', 'qstrike', 'qdip', 'qrake']:
                if self.ehead[key] != -99.0:
                    evdata[key] = np.round(self.ehead[key], 6)
            for imag in range(1, 4):
                key1, key2 = 'qmag{:}'.format(imag), 'qmag{:}type'.format(imag)
                if (self.ehead[key1] != 0.0) or (self.ehead[key2] != '    '):
                    evdata[key1] = np.round(self.ehead[key1], 3)
                    evdata[key2] = self.ehead[key2].strip()

        # Loop over time
        for i, wf in enumerate(self.waveforms):

            # Initialize stats for trace
            stats = {}

            # Assemble mandatory header information
            stats['delta'] = wf['dt']
            stats['sampling_rate'] = 1 / wf['dt']
            if wf['gain'] > 0:
                stats['calib'] = wf['gain']
            else:
                stats['calib'] = 1.0
            stats['npts'] = wf['npts']
            stats['network'] = wf['stype'].strip()
            stats['station'] = wf['stname'].strip()
            stats['location'] = wf['loccode'].strip()
            if stats['location'] == '--':
                stats['location'] = ''
            stats['channel'] = wf['chnm'].strip()
            stats['starttime'] = UTCDateTime(wf['syr'], wf['smon'], wf['sdy'],
                                             wf['shr'], wf['smn'], np.round(wf['ssc'], 5))
            
            # IV DEBUG
            # ID = '.'.join([stats['network'],stats['station'],stats['location'],stats['channel']])

            # Optional: preserve event information in header
            if keep_evdata:
                stats['event_data'] = evdata

            # Optional: preserve station information in header
            if keep_stdata:
                stats['station_data'] = {}
                for key in ['compazi', 'compang', 'deldist', 'sazi', 'qazi', 'slat', 'slon', 'selev']:
                    if wf[key] != -99.0:
                        stats['station_data'][key] = np.round(wf[key], 6)

            # Optional: preserve pick information in header
            if keep_pkdata:
                stats['pick_data'] = {}
                stats['pick_data']['tdif'] = np.round(wf['tdif'], 6)
                stats['pick_data']['ppolarity'] = wf['ppolarity'].strip()
                for ipick in range(1, 5):
                    key1 = 'pick{:}'.format(ipick)
                    
                    # IV NOTE: causes issues if trace starts after pick time. changed from 'if wf[key1] > 0:'
                    if wf[key1] != 0 and wf[key1] != -99.0:
                        key2, key3 = 'pick{:}name'.format(ipick), 'pick{:}q'.format(ipick)
                        stats['pick_data'][key1] = np.round(wf[key1], 6)
                        stats['pick_data'][key2] = wf[key2].strip()
                        stats['pick_data'][key3] = wf[key3].strip()
                # if ID=='NP.1809..HNE': print(stats['pick_data'])

            # Update stream
            st += Trace(data = wf['data'], header = stats)
            
        # return Obspy stream
        return st

        ###########################################################

    ### Function to Convert From ObsPy Stream to EFS
    # (todo: bytepos?)
    def from_obspy(st, evhead={}, invs={}):

        '''
        Function to return EFS from obspy stream.
        Optional argument for event header to populate EFS event header field.
        '''

        # Error checking
        try:
            from obspy.core import UTCDateTime
            from obspy import geodetics
        except:
            print('EFS.from_obspy() is not available. ObsPy is not installed.')
            raise

        # initialize EFS: blank fhead, ehead, waveforms
        #    choose data type from first trace
        if (st[0].data[0].dtype=="float64")or(st[0].data[0].dtype=="int64"):
            efs_data = EFS(None,prec_wf=64)
        else:
            efs_data = EFS(None,prec_wf=32)

        # set default tshead
        efs_data.fhead['bytetype'] = 1
        efs_data.fhead['eheadtype'] = 1
        efs_data.fhead['nbytes_ehead'] = 264
        efs_data.fhead['tsheadtype'] = 1
        efs_data.fhead['nbytes_tshead'] = 268

        # set up ehead, reading from evhead where possible
        for key in ['efslabel', 'datasource']:
            if key in evhead:
                efs_data.ehead[key] = '{:40s}'.format(evhead[key])
            else:
                efs_data.ehead[key] = '                                        '
        for key in ['qtype', 'qmag1type', 'qmag2type', 'qmag3type',
                    'qmomenttype', 'qlocqual', 'qfocalqual']:
            if key in evhead:
                efs_data.ehead[key] = '{:4s}'.format(evhead[key])
            else:
                efs_data.ehead[key] = '    '
        for key in ['qlat', 'qlon', 'qdep', 'qsc', 'qmag1', 'qmag2', 'qmag3',
                    'qmoment', 'qstrike', 'qdip', 'qrake', 'qyr', 'qmon',
                    'qdy', 'qhr', 'qmn', 'cuspid']:
            if key in evhead:
                efs_data.ehead[key] = evhead[key]
            elif key in ['qlat', 'qlon', 'qdep', 'qsec', 'qmoment']:
                efs_data.ehead[key] = -999.
            else:
                efs_data.ehead[key] = -999

        # loop over traces in stream
        efs_data.ehead['numts'] = len(st)
        efs_data.ehead['maxnumts'] = len(st)
        for ii, tr in enumerate(st):

            # tshead: blank fields
            wf = {
                'datasource': '        ', 'sensor': '        ', 'units': '        ', 'dva': '    ',
                'pick1q': '    ', 'pick2q': '    ', 'pick3q': '    ', 'pick4q': '    ',
                'pick1name': '    ', 'pick2name': '    ', 'pick3name': '    ', 'pick4name': '    ',
                'ppolarity': '    ', 'compazi': -99.0, 'compang': -99.0, 'f1': -1.0, 'f2': -1.0,
                'deldist': -99.0, 'slat': -99.0, 'slon': -99.0, 'selev': -99.0,
                'sazi': -99.0, 'qazi': -99.0, 'pick1': 0.0, 'pick2': 0.0, 'pick3': 0.0, 'pick4': 0.0,
                'problem': '    '
            }

            # tshead: fields from stats
            # only station lat, lon, and elevation taken from inventory
            if 'invs' in locals():
                # print(len(invs))        #IV edit
                tmp1 = invs.select(station = tr.stats.station)
                wf['slat'] = tmp1[0].stations[0].latitude
                wf['slon'] = tmp1[0].stations[0].longitude
                wf['selev'] = tmp1[0].stations[0].elevation

                tmp2 = geodetics.gps2dist_azimuth(wf['slat'], wf['slon'], evhead['qlat'], evhead['qlon'])
                wf['qazi'] = tmp2[2]
                wf['sazi'] = tmp2[1]
                wf['deldist'] = tmp2[0] / 1E3
                # wf['deldist'] = geodetics.kilometer2degrees(tmp2[0] / 1E3) # IV replaced with above line 22 Nov 2021- EFS field should be in km in documentation

            wf['dt'] = tr.stats.delta
            wf['syr'] = tr.stats.starttime.year
            wf['smon'] = tr.stats.starttime.month
            wf['sdy'] = tr.stats.starttime.day
            wf['shr'] = tr.stats.starttime.hour
            wf['smn'] = tr.stats.starttime.minute
            wf['ssc'] = float(tr.stats.starttime.second + tr.stats.starttime.microsecond / 1.e6)
            wf['gain'] = tr.stats.calib
            wf['npts'] = tr.stats.npts

            # net, sta, chan, loc
            wf['stype'] = "{:4s}".format(tr.stats.network)
            wf['stname'] = "{:8s}".format(tr.stats.station)
            wf['loccode'] = "{:8s}".format(tr.stats.location)
            wf['chnm'] = "{:4s}".format(tr.stats.channel)

            # calculate tdif = origin time - starttime
            try:
                qtime = UTCDateTime(
                    efs_data.ehead['qyr'], efs_data.ehead['qmon'], efs_data.ehead['qdy'],
                    efs_data.ehead['qhr'], efs_data.ehead['qmn'], efs_data.ehead['qsc']
                )
                wf['tdif'] = tr.stats.starttime - qtime # IV switched order to switch sign
            except:
                wf['tdif'] = 0.0

            # waveform data
            wf['data'] = tr.data

            # add to list
            efs_data.waveforms.append(wf)

        # return
        return efs_data


### Function to write EFS file
# file is EFSPATH + efsname
# data comes from efs_data
# prec_wf is the time series array precision (default is f32)
# prec_bp is the byte position array precision (default is i32)
def export_efs(EFSPATH, efsname, efs_data, prec_wf = np.float32, prec_bp = "i"):
    
    # open file
    fname2 = EFSPATH + efsname
    f22 = open(fname2, "wb")

    # write file header (20 bytes)
    f22.write(struct.pack("i", efs_data.fhead['bytetype']))
    f22.write(struct.pack("i", efs_data.fhead['eheadtype']))
    f22.write(struct.pack("i", efs_data.fhead['nbytes_ehead']))
    f22.write(struct.pack("i", efs_data.fhead['tsheadtype']))
    f22.write(struct.pack("i", efs_data.fhead['nbytes_tshead']))

    # write event header (184 bytes)
    f22.write(struct.pack("40s", efs_data.ehead['efslabel'].encode()))
    f22.write(struct.pack("40s", efs_data.ehead['datasource'].encode()))
    f22.write(struct.pack("i", efs_data.ehead['maxnumts']))
    f22.write(struct.pack("i", efs_data.ehead['numts']))
    f22.write(struct.pack("i", efs_data.ehead['cuspid']))
    f22.write(struct.pack("4s", efs_data.ehead['qtype'].encode()))
    f22.write(struct.pack("4s", efs_data.ehead['qmag1type'].encode()))
    f22.write(struct.pack("4s", efs_data.ehead['qmag2type'].encode()))
    f22.write(struct.pack("4s", efs_data.ehead['qmag3type'].encode()))
    f22.write(struct.pack("4s", efs_data.ehead['qmomenttype'].encode()))
    f22.write(struct.pack("4s", efs_data.ehead['qlocqual'].encode()))
    f22.write(struct.pack("4s", efs_data.ehead['qfocalqual'].encode()))
    f22.write(struct.pack("f", efs_data.ehead['qlat']))
    f22.write(struct.pack("f", efs_data.ehead['qlon']))
    f22.write(struct.pack("f", efs_data.ehead['qdep']))
    f22.write(struct.pack("f", efs_data.ehead['qsc']))
    f22.write(struct.pack("f", efs_data.ehead['qmag1']))
    f22.write(struct.pack("f", efs_data.ehead['qmag2']))
    f22.write(struct.pack("f", efs_data.ehead['qmag3']))
    f22.write(struct.pack("f", efs_data.ehead['qmoment']))
    f22.write(struct.pack("f", efs_data.ehead['qstrike']))
    f22.write(struct.pack("f", efs_data.ehead['qdip']))
    f22.write(struct.pack("f", efs_data.ehead['qrake']))
    f22.write(struct.pack("i", efs_data.ehead['qyr']))
    f22.write(struct.pack("i", efs_data.ehead['qmon']))
    f22.write(struct.pack("i", efs_data.ehead['qdy']))
    f22.write(struct.pack("i", efs_data.ehead['qhr']))
    f22.write(struct.pack("i", efs_data.ehead['qmn']))

    # 20 4-byte field (80 bytes)
    for idum in range(0, 20):
        f22.write(struct.pack("i", idum * 0))

    # byte positions for all time series
    for ipos in range(0, efs_data.ehead['numts']):
        f22.write(struct.pack(prec_bp, efs_data.ehead['bytepos'][ipos]))

    # ilen = 100000
    # np.array(np.zeros(ilen*efs_data.ehead['numts']), dtype=np.uint32).tofile(f22)

    # write time series data
    for ii in range(0, efs_data.ehead['numts']):
        f22.seek(efs_data.ehead['bytepos'][ii])
        f22.write(struct.pack("8s", efs_data.waveforms[ii]['stname'].encode()))
        f22.write(struct.pack("8s", efs_data.waveforms[ii]['loccode'].encode()))
        f22.write(struct.pack("8s", efs_data.waveforms[ii]['datasource'].encode()))
        f22.write(struct.pack("8s", efs_data.waveforms[ii]['sensor'].encode()))
        f22.write(struct.pack("8s", efs_data.waveforms[ii]['units'].encode()))
        f22.write(struct.pack("4s", efs_data.waveforms[ii]['chnm'].encode()))
        f22.write(struct.pack("4s", efs_data.waveforms[ii]['stype'].encode()))
        f22.write(struct.pack("4s", efs_data.waveforms[ii]['dva'].encode()))
        f22.write(struct.pack("4s", efs_data.waveforms[ii]['pick1q'].encode()))
        f22.write(struct.pack("4s", efs_data.waveforms[ii]['pick2q'].encode()))
        f22.write(struct.pack("4s", efs_data.waveforms[ii]['pick3q'].encode()))
        f22.write(struct.pack("4s", efs_data.waveforms[ii]['pick4q'].encode()))
        f22.write(struct.pack("4s", efs_data.waveforms[ii]['pick1name'].encode()))
        f22.write(struct.pack("4s", efs_data.waveforms[ii]['pick2name'].encode()))
        f22.write(struct.pack("4s", efs_data.waveforms[ii]['pick3name'].encode()))
        f22.write(struct.pack("4s", efs_data.waveforms[ii]['pick4name'].encode()))
        f22.write(struct.pack("4s", efs_data.waveforms[ii]['ppolarity'].encode()))
        f22.write(struct.pack("4s", efs_data.waveforms[ii]['problem'].encode()))
        f22.write(struct.pack("i", efs_data.waveforms[ii]['npts']))
        f22.write(struct.pack("i", efs_data.waveforms[ii]['syr']))
        f22.write(struct.pack("i", efs_data.waveforms[ii]['smon']))
        f22.write(struct.pack("i", efs_data.waveforms[ii]['sdy']))
        f22.write(struct.pack("i", efs_data.waveforms[ii]['shr']))
        f22.write(struct.pack("i", efs_data.waveforms[ii]['smn']))
        f22.write(struct.pack("f", efs_data.waveforms[ii]['compazi']))
        f22.write(struct.pack("f", efs_data.waveforms[ii]['compang']))
        f22.write(struct.pack("f", efs_data.waveforms[ii]['gain']))
        f22.write(struct.pack("f", efs_data.waveforms[ii]['f1']))
        f22.write(struct.pack("f", efs_data.waveforms[ii]['f2']))
        f22.write(struct.pack("f", efs_data.waveforms[ii]['dt']))
        f22.write(struct.pack("f", efs_data.waveforms[ii]['ssc']))
        f22.write(struct.pack("f", efs_data.waveforms[ii]['tdif']))
        f22.write(struct.pack("f", efs_data.waveforms[ii]['slat']))
        f22.write(struct.pack("f", efs_data.waveforms[ii]['slon']))
        f22.write(struct.pack("f", efs_data.waveforms[ii]['selev']))
        f22.write(struct.pack("f", efs_data.waveforms[ii]['deldist']))
        f22.write(struct.pack("f", efs_data.waveforms[ii]['sazi']))
        f22.write(struct.pack("f", efs_data.waveforms[ii]['qazi']))
        f22.write(struct.pack("f", efs_data.waveforms[ii]['pick1']))
        f22.write(struct.pack("f", efs_data.waveforms[ii]['pick2']))
        f22.write(struct.pack("f", efs_data.waveforms[ii]['pick3']))
        f22.write(struct.pack("f", efs_data.waveforms[ii]['pick4'] * 0 - 99))
        
        # dummy spots
        for idum in range(0, 20):
            f22.write(struct.pack("i", idum * 0))

        # waveform arrays
        tmp1 = efs_data.waveforms[ii]['data']
        #np.array(tmp1, dtype = itstype).tofile(f22)
        tmp1.astype(prec_wf).tofile(f22)

    f22.close() # close

### Simple function to convert ObsPy catalog to EFS event header
def cat2ehead(st1,cat):
    import obspy
    ehead = {}
    ehead['efslabel'] = "{:<40}".format(' ')
    ehead['datasource'] = "{:<40}".format('miniSEED')
    ehead['maxnumts'] = 1000
    ehead['numts'] = len(st1)
    ehead['cuspid'] = 0
    ehead['qtype'] = "{:<4}".format(' ')
    ehead['qmag1type'] = cat[0].magnitudes[0].magnitude_type
    ehead['qmag2type'] = "{:<4}".format(' ')
    ehead['qmag3type'] = "{:<4}".format(' ')
    ehead['qmomenttype'] = "{:<4}".format(' ')
    ehead['qlocqual'] = "{:<4}".format(' ')
    ehead['qfocalqual'] = "{:<4}".format(' ')
    ehead['qlat'] = cat[0].origins[0].latitude
    ehead['qlon'] = cat[0].origins[0].longitude
    ehead['qdep'] = cat[0].origins[0].depth
    ehead['qsc'] = cat[0].origins[0].time.second
    ehead['qmag1'] = cat[0].magnitudes[0].mag
    ehead['qmag2'] = 0
    ehead['qmag3'] = 0
    ehead['qmoment'] = 0
    ehead['qstrike'] = 0
    ehead['qdip'] = 0
    ehead['qrake'] = 0
    ehead['qyr'] = cat[0].origins[0].time.year
    ehead['qmon'] = cat[0].origins[0].time.month
    ehead['qdy'] = cat[0].origins[0].time.day
    ehead['qhr'] = cat[0].origins[0].time.hour
    ehead['qmn'] = cat[0].origins[0].time.minute

    return ehead

# IV edits below

#### EFS_head Class Definition ###
class EFS_head():
    '''
    Class definition for EFS-format data.
    Basic initialization syntax: edata = EFS("path_to_efs_file").
    '''

    ### --- Initialization --- ###
    #  efsfname: name of file to read
    #  prec_wf: precision for waveform arrays, 32 or 64
    #  prec_bp: precision for byteposition arrays
    def __init__(self, efsfname=None, prec_wf=np.float32, prec_bp=np.int32):

        # initialize fields
        self.fhead = {}
        self.ehead = {}
        self.waveforms = []

        # return here without file
        if efsfname is None:
            return

        # Open the EFS binary file
        f = open(efsfname, 'rb')

        # Assemble file header
        self.fhead['bytetype'] = struct.unpack('i', f.read(4))[0]
        self.fhead['eheadtype'] = struct.unpack('i', f.read(4))[0]
        self.fhead['nbytes_ehead'] = struct.unpack('i', f.read(4))[0]
        self.fhead['tsheadtype'] = struct.unpack('i', f.read(4))[0]
        self.fhead['nbytes_tshead'] = struct.unpack('i', f.read(4))[0]

        # Assemble event header
        self.ehead['efslabel'] = f.read(40).decode('UTF-8')
        self.ehead['datasource'] = f.read(40).decode('UTF-8')
        self.ehead['maxnumts'] = struct.unpack('i', f.read(4))[0]
        self.ehead['numts'] = struct.unpack('i', f.read(4))[0]
        self.ehead['cuspid'] = struct.unpack('i', f.read(4))[0]
        self.ehead['qtype'] = f.read(4).decode('UTF-8')
        self.ehead['qmag1type'] = f.read(4).decode('UTF-8')
        self.ehead['qmag2type'] = f.read(4).decode('UTF-8')
        self.ehead['qmag3type'] = f.read(4).decode('UTF-8')
        self.ehead['qmomenttype'] = f.read(4).decode('UTF-8')
        self.ehead['qlocqual'] = f.read(4).decode('UTF-8')
        self.ehead['qfocalqual'] = f.read(4).decode('UTF-8')
        self.ehead['qlat'] = struct.unpack('f', f.read(4))[0]
        self.ehead['qlon'] = struct.unpack('f', f.read(4))[0]
        self.ehead['qdep'] = struct.unpack('f', f.read(4))[0]
        self.ehead['qsc'] = struct.unpack('f', f.read(4))[0]
        self.ehead['qmag1'] = struct.unpack('f', f.read(4))[0]
        self.ehead['qmag2'] = struct.unpack('f', f.read(4))[0]
        self.ehead['qmag3'] = struct.unpack('f', f.read(4))[0]
        self.ehead['qmoment'] = struct.unpack('f', f.read(4))[0]
        self.ehead['qstrike'] = struct.unpack('f', f.read(4))[0]
        self.ehead['qdip'] = struct.unpack('f', f.read(4))[0]
        self.ehead['qrake'] = struct.unpack('f', f.read(4))[0]
        self.ehead['qyr'] = struct.unpack('i', f.read(4))[0]
        self.ehead['qmon'] = struct.unpack('i', f.read(4))[0]
        self.ehead['qdy'] = struct.unpack('i', f.read(4))[0]
        self.ehead['qhr'] = struct.unpack('i', f.read(4))[0]
        self.ehead['qmn'] = struct.unpack('i', f.read(4))[0]

        # 20 4-byte fields reserved for future uses - skip
        for idum in range(0, 20):
            dummy = struct.unpack('i', f.read(4))[0]

        # Get byte positions for all time series (i64 or i32 arrays)
        bytepos = np.fromfile(f, dtype = prec_bp, count = self.ehead['numts'])
        self.ehead['bytepos'] = bytepos

        # Now loop over all the time series
        for ii in range(0, len(bytepos)):

            # Assemble tshead
            f.seek(bytepos[ii])
            tshead = {}
            tshead['stname'] = f.read(8).decode('UTF-8')
            tshead['loccode'] = f.read(8).decode('UTF-8')
            tshead['datasource'] = f.read(8).decode('UTF-8')
            tshead['sensor'] = f.read(8).decode('UTF-8')
            tshead['units'] = f.read(8).decode('UTF-8')
            tshead['chnm'] = f.read(4).decode('UTF-8')
            tshead['stype'] = f.read(4).decode('UTF-8')
            tshead['dva'] = f.read(4).decode('UTF-8')
            tshead['pick1q'] = f.read(4).decode('UTF-8')
            tshead['pick2q'] = f.read(4).decode('UTF-8')
            tshead['pick3q'] = f.read(4).decode('UTF-8')
            tshead['pick4q'] = f.read(4).decode('UTF-8')
            tshead['pick1name'] = f.read(4).decode('UTF-8')
            tshead['pick2name'] = f.read(4).decode('UTF-8')
            tshead['pick3name'] = f.read(4).decode('UTF-8')
            tshead['pick4name'] = f.read(4).decode('UTF-8')
            tshead['ppolarity'] = f.read(4).decode('UTF-8')
            tshead['problem'] = f.read(4).decode('UTF-8')
            tshead['npts'] = struct.unpack('i', f.read(4))[0]
            tshead['syr'] = struct.unpack('i', f.read(4))[0]
            tshead['smon'] = struct.unpack('i', f.read(4))[0]
            tshead['sdy'] = struct.unpack('i', f.read(4))[0]
            tshead['shr'] = struct.unpack('i', f.read(4))[0]
            tshead['smn'] = struct.unpack('i', f.read(4))[0]
            tshead['compazi'] = struct.unpack('f', f.read(4))[0]
            tshead['compang'] = struct.unpack('f', f.read(4))[0]
            tshead['gain'] = struct.unpack('f', f.read(4))[0]
            tshead['f1'] = struct.unpack('f', f.read(4))[0]
            tshead['f2'] = struct.unpack('f', f.read(4))[0]
            tshead['dt'] = struct.unpack('f', f.read(4))[0]
            tshead['ssc'] = struct.unpack('f', f.read(4))[0]
            tshead['tdif'] = struct.unpack('f', f.read(4))[0]
            tshead['slat'] = struct.unpack('f', f.read(4))[0]
            tshead['slon'] = struct.unpack('f', f.read(4))[0]
            tshead['selev'] = struct.unpack('f', f.read(4))[0]
            tshead['deldist'] = struct.unpack('f', f.read(4))[0]
            tshead['sazi'] = struct.unpack('f', f.read(4))[0]
            tshead['qazi'] = struct.unpack('f', f.read(4))[0]
            tshead['pick1'] = struct.unpack('f', f.read(4))[0]
            tshead['pick2'] = struct.unpack('f', f.read(4))[0]
            tshead['pick3'] = struct.unpack('f', f.read(4))[0]
            tshead['pick4'] = struct.unpack('f', f.read(4))[0]

            # 20 4-byte fields reserved for future uses - skip
            # for idum in range(0, 20):
            #     dummy = struct.unpack('i', f.read(4))[0]

            # Read the time-series itself
            # data = np.fromfile(f, dtype = prec_wf, count = tshead['npts'])

            # Bundle tsheader and time-series for this waveform into efsdata, then append to list
            efsdata = tshead
            # efsdata['data'] = data
            self.waveforms.append(efsdata)