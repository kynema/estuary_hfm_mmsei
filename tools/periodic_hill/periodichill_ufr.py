import numpy as np
import subprocess, sys, os
import argparse
import ruamel.yaml as yaml
import time
import netCDF4 as nc
import matplotlib.pyplot as plt
import pandas as pd

def plot_probes(data,name):

    time = data.variables['time'][:]
    vel_x = data['probe1']['velocityx']
    vel_y = data['probe1']['velocityy']
    vel_z = data['probe1']['velocityz']
    nprobes = data['probe1'].dimensions['num_points'].size

    fig, ax = plt.subplots(3, 1, figsize=(12, 8), sharex="col")

    for i in range(nprobes):
        ax[0].plot(time,vel_x[:,i], label="Point "+str(i))
        ax[0].set_title("Velocity Timeseries")
        ax[0].set_ylabel("X Velocity")
        ax[0].grid(True, which="both", linestyle="--", linewidth=0.1)

        ax[1].plot(time,vel_y[:,i], label="Point "+str(i))
        ax[1].set_ylabel("Y Velocity")
        ax[1].grid(True, which="both", linestyle="--", linewidth=0.1)

        ax[2].plot(time,vel_y[:,i], label="Point "+str(i))
        ax[2].set_ylabel("Z Velocity")
        ax[2].set_xlabel("Time")
        ax[2].grid(True, which="both", linestyle="--", linewidth=0.1)
        ax[2].legend()

    plt.tight_layout()
    plt.savefig('probe_sampler_timeseries.png')

def load_cases(cases):
    planedata = {}
    pointdata = {}
    

    for c in cases:
        postproc_dir = os.path.join(c['path'],'post_processing')
        planedata[c['name']] = nc.Dataset(os.path.join(postproc_dir, c['plane_sample_file']))
        pointdata[c['name']] = nc.Dataset(os.path.join(postproc_dir, c['point_sample_file']))

    return planedata,pointdata

def load_exp(exp):
    
    pardir = exp['parent_dir']
    variables = exp['variables']
    profile = pd.read_csv(exp['profile'])

    alldata = {}

    for e in exp['sources']:
        x = e['xval']
        alldata[x] = {}
        filepath = os.path.join(pardir,e['files'][0])
        # Variables: y/h; u/u_b; v/u_b; u'u'/u_b^2; v'v'/u_b^2; u'v'/u_b^2
        alldata[x] = pd.read_csv(filepath, sep=',', skiprows=6, header=None)

    return alldata,profile

def make_means(planedata,expdata,profile,name,opt,label):

    vel_x = planedata['plane1']['velocityx']
    vel_y = planedata['plane1']['velocityy']
    vel_z = planedata['plane1']['velocityz']

    #print(planedata['plane1'])
    print("Len of vel_x:",len(vel_x))

    nsteps = planedata['time'].shape[0]
    nwidth = planedata['plane1'].ijk_dims[0]
    nlength = planedata['plane1'].ijk_dims[1]
    noffsets = planedata['plane1'].ijk_dims[2]
    ntp = nwidth*nlength
    offsets = planedata['plane1'].offsets

    #print(name,"nsteps",nsteps)

    ms = opt[0]
    me = opt[1] #nsteps

    #print(planedata['plane1'].axis1[1])
    #print(planedata['plane1'].axis2[2])


    yvals = np.linspace(0,planedata['plane1'].axis1[1], num=nwidth, endpoint=True)
    zvals = np.linspace(0,planedata['plane1'].axis2[2], num=nlength, endpoint=True)
    Z,Y = np.meshgrid(zvals,yvals)

    nmean = len(range(ms,me))
    thismean = np.zeros((nwidth,nlength,noffsets))

    # Calculate temporal mean
    for step in range(ms,me):
        v = vel_x[step,:].reshape((nwidth,nlength,noffsets),order='F')
        thismean = thismean + v/nmean

    # Calculate spanwise
    hmean = []
    for o in range(len(offsets)):
        hmean.append(thismean[:,:,o].reshape(nwidth,nlength).mean(axis=0))

    centerslice = []
    for o in range(len(offsets)):
        centerslice.append(thismean[:,:,o].reshape(nwidth,nlength)[(nwidth // 2),:])

    maxumid = np.max(hmean[5][:])

    return {'t_mean': thismean,
            'h_mean': hmean,
            'center_slice':centerslice,
            'Zmesh': Z,
            'Ymesh': Y,
            'zvals': zvals,
            'ns':nsteps,
            'nw':nwidth,
            'nl':nlength,
            'noffsets':noffsets,
            'nmean':nmean,
            'offsets':offsets,
            'vel_x': vel_x,
            'vel_y': vel_y,
            'vel_z': vel_z,
            'opt':opt,
            'label': label,
            'ummid': maxumid
            }


def plot_planes(m,expdata,profile):

    for name in m:
        # Plot Spanwise Mean
        print('Plotting:',name)
        fig, ax = plt.subplots(1, 1, figsize=(12, 8), sharex="col")
        ax.plot(profile['x'],profile['z']-0.00132812,label="Hill Geometry",color="red")
        ax.fill_between(profile['x'],profile['z']-0.00132812, color="red", alpha=0.7)
        scale = m[name]['opt'][3]
        u_b = m[name]['opt'][4]

        for o in range(m[name]['noffsets']):
            ax.plot(m[name]['offsets'][o]+m[name]['h_mean'][o]/u_b/scale,m[name]['zvals'],label="Offset "+str(o))

        for i,e in enumerate(expdata):
            z_h = expdata[e]['U'][0]*0.028
            U_m = expdata[e]['U'][1]/scale
            ax.plot(U_m+float(e),z_h,label="Exp",color="black")
        
        plt.axis('equal')
        plt.title("Mean Velocity")
        plt.xlabel("Velocity Profile")
        plt.ylabel("Z axis")
        plt.savefig('plane_u_mean_'+name+'.png')


def plot_planes_all(m,expdata,profile):

    h = 0.028

    uls = ['-',':','--','-.',(0, (5, 5))]

    cp = [
        "#E69F00",  # Orange
        "#56B4E9",  # Sky Blue
        "#009E73",  # Bluish Green
        "#F0E442",  # Yellow
        "#0072B2",  # Blue
        "#D55E00",  # Vermillion
        "#CC79A7",  # Reddish Purple
        "#999999",  # Gray
        "#E41A1C",  # Red (colorblind safe variant)
        "#377EB8",  # Blue variant
        "#A6761D"   # Brownish (added for 11th color)
    ]
    
    ptype = ['h_mean','center_slice']
    
    for pt in ptype:
        fig, ax = plt.subplots(1, 1, figsize=(12, 8), sharex="col")
        scale = 1.0
        
        for l,name in enumerate(m):
            print('Max U:',m[name]['ummid'])
            u_b = m[name]['ummid']*0.85
            for o in range(m[name]['noffsets']):
                ax.plot(m[name]['offsets'][o]/h+m[name][pt][o]/scale/u_b,m[name]['zvals']/h,label=m[name]['label'],linestyle=uls[l+1],color=cp[o])

        # Variables: y/h; u/u_b; v/u_b; u'u'/u_b^2; v'v'/u_b^2; u'v'/u_b^2
        for i,e in enumerate(expdata):
            z_h = expdata[e][0]
            U_m = expdata[e][1]/scale
            ax.plot(U_m+float(e)/h,z_h,label="Exp",color="black",linestyle=uls[0])

        ax.plot(profile['x']/h,profile['z']/h,label="Hill Geometry",color="#222222")
        ax.fill_between(profile['x']/h,profile['z']/h, color="#222222", alpha=0.8)
        
        handles, labels = ax.get_legend_handles_labels()

        grouped_handles = {}
        grouped_labels = {}

        for handle, label in zip(handles, labels):
            ls = handle.get_linestyle()
            if (ls not in grouped_handles):
                grouped_handles[ls] = handle
                grouped_labels[ls] = label # Keeps the first label encountered for this style

        ax.legend(handles=grouped_handles.values(), labels=grouped_labels.values())

        plt.axis('equal')
        plt.title("Scaled Mean Velocity")
        plt.xlabel("x/h with Mean Velocity (U/U_b)/100")
        plt.ylabel("Z (m)")
        plt.savefig('ufr_' + pt + '_flow.png')
    

def main():

    parser = argparse.ArgumentParser(description="Analyze Periodic Hill")

    parser.add_argument(
        "-i",
        "--infile",
        help="Input YAML file (must be present in the current directory)",
        required=True,
        type=str,
    )

    args = parser.parse_args()

    # Load input yaml
    yml = yaml.YAML(typ='safe', pure=True)    
    with open(args.infile, 'r') as stream:
        inp = yml.load(stream)


    ###############################################################
    # YAML inputs:
    
    cases = inp['cases']
    exp = inp['data']

    ###############################################################


    planedata,pointdata = load_cases(cases)

    expdata,profile = load_exp(exp)

    output = {}

    plt.rcParams['font.size'] = 16
        
    for c in cases:
        print("Processing:",c['name'])
        if c['probe_plots']:
            plot_probes(pointdata[c['name']],c['name'])
        if c['plane_processing']:
            output[c['name']] = make_means(planedata[c['name']],expdata,profile,c['name'],c['opt'],c['label'])
            print("Plane output steps:",output[c['name']]['ns'])
    
    plot_planes_all(output,expdata,profile)

        

if __name__=="__main__":
    main()