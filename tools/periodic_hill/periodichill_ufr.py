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
    zcoords = data['probe1'].variables['coordinates']

    fig, ax = plt.subplots(3, 1, figsize=(12, 8), sharex="col")

    for i in range(nprobes):
        #print("Probe z=",zcoords[i][2],"Mean X Velocity:",np.mean(vel_x[:,i])) 

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

def load_ugf_planes(sample_paths):
    # Load UGF plane data from multiple files and combine into a single dataset
    combined_data = []
    for path in sample_paths:
        data = pd.read_csv(path, sep=r'\s+', skipinitialspace=True,skiprows=1)
        combined_data.append(data)

    return combined_data

def load_cases(cases):
    planedata = {}
    pointdata = {}
    
    for c in cases:
        if c['type'] == 'sgf':
            postproc_dir = os.path.join(c['path'],'post_processing')
            planedata[c['name']] = nc.Dataset(os.path.join(postproc_dir, c['plane_sample_file'][0]))
            pointdata[c['name']] = nc.Dataset(os.path.join(postproc_dir, c['point_sample_file']))
        else:
            postproc_dir = os.path.join(c['path'])
            popt = c['plane_options']
            plane_i = range(popt[0],popt[1]+popt[2],popt[2])
            sample_paths = [os.path.join(postproc_dir, c['plane_template'].replace('{i}',f"{sp:07}")) for sp in plane_i]
            planedata[c['name']] = load_ugf_planes(sample_paths) # List of flat files for each time output
            pointdata[c['name']] = []

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

def make_means_sgf(planedata,name,opt,label):

    vel_x = planedata['plane1']['velocityx']
    vel_y = planedata['plane1']['velocityy']
    vel_z = planedata['plane1']['velocityz']

    #print(planedata['plane1'])
    print("Total flowthroughs in vel_x:",len(vel_x))

    nsteps = planedata['time'].shape[0]
    nwidth = planedata['plane1'].ijk_dims[0]
    nlength = planedata['plane1'].ijk_dims[1]
    noffsets = planedata['plane1'].ijk_dims[2]
    ntp = nwidth*nlength
    offsets = planedata['plane1'].offsets

    ms = opt[0]
    me = opt[1] #nsteps

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

    maxumid = np.max(hmean[0][:])

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

def zero_under_hill(hilldata,meandata):

    # Zero out the mean data under the hill geometry
    hill_z = hilldata['z']
    hill_x = hilldata['x']

    for i in range(len(meandata)):
        x_val = meandata.iloc[i]['coordinates[0]']
        z_val = meandata.iloc[i]['coordinates[2]']

        # Find the corresponding hill height at this x position
        hill_height = np.interp(x_val, hill_x, hill_z)

        if z_val < hill_height:
            meandata.at[i, 'velocity_probe[0]'] = 0.0
            meandata.at[i, 'velocity_probe[1]'] = 0.0
            meandata.at[i, 'velocity_probe[2]'] = 0.0

    return meandata

def make_means_ugf(planedata,name,opt,label):

    time_mean = sum(planedata) / len(planedata)

    hilldata = pd.read_csv('periodic_hill_ufr.csv')
    time_mean = zero_under_hill(hilldata,time_mean)

    # time_mean indices
    # noffsets is slowest changing, then Index J, then Index I is fastest changing
    # j=length, i=width, o=offsets

    nsteps = len(planedata)
    nwidth = np.unique(planedata[0]['Index_i']).size
    nlength = np.unique(planedata[0]['Index_j']).size
    noffsets = np.unique(planedata[0]['#Plane_Number']).size
    ntp = nwidth*nlength
    offsets = np.unique(planedata[0]['coordinates[0]'])
    nmean = len(time_mean)

    vel_x = []
    vel_y = []
    vel_z = []

    for data in planedata:
        vel_x.append(np.array(data['velocity_probe[0]']))
        vel_y.append(np.array(data['velocity_probe[1]']))
        vel_z.append(np.array(data['velocity_probe[2]']))

    ms = opt[0]
    me = opt[1] #nsteps

    yvals = np.unique(planedata[0]['coordinates[1]'])
    zvals = np.unique(planedata[0]['coordinates[2]'])
    Z,Y = np.meshgrid(zvals,yvals)

    # Calculate temporal mean
    thismean = np.array(time_mean['velocity_probe[0]']).reshape(noffsets,nlength,nwidth, order='C').transpose(2,1,0) # Reshape to (noffsets,nwidth,nlength)

    #print(np.array(thismean))

    # Calculate spanwise
    hmean = []
    for o in range(len(offsets)):
        hmean.append(thismean[:,:,o].reshape(nwidth,nlength).mean(axis=0))

    centerslice = []
    for o in range(len(offsets)):
        centerslice.append(thismean[:,:,o].reshape(nwidth,nlength)[(nwidth // 2),:])

    maxumid = np.max(hmean[0][:])

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

    uls = ['--','-.',':','-','-']

    cp = [
        "#D55E00",  # Vermillion
        "#E69F00",  # Orange
        "#CC79A7",  # Reddish Purple
        "#009E73",  # Bluish Green
        "#56B4E9",  # Sky Blue
        "#F0E442",  # Yellow
        "#0072B2",  # Blue
        "#999999",  # Gray
        "#E41A1C",  # Red (colorblind safe variant)
        "#377EB8",  # Blue variant
        "#A6761D"   # Brownish (added for 11th color)
    ]
    
    ptype = ['h_mean','center_slice']
    
    for pt in ptype:
        fig, ax = plt.subplots(1, 1, figsize=(14, 8), sharex="col")
        scale = 1.0
        
        for l,name in enumerate(m):
            u_b = m[name]['ummid']/1.0593
            print("U_bulk for",name,"=",u_b)
            for o in range(m[name]['noffsets']):
                ax.plot(m[name]['offsets'][o]/h+m[name][pt][o]/scale/u_b,m[name]['zvals']/h,label=m[name]['label'],linestyle=uls[l],color=cp[l],lw=1.5)

        # Variables: y/h; u/u_b; v/u_b; u'u'/u_b^2; v'v'/u_b^2; u'v'/u_b^2
        for i,e in enumerate(expdata):
            z_h = expdata[e][0]
            U_m = expdata[e][1]/scale
            ax.plot(U_m+float(e)/h,z_h,label="Exp",color="black",linestyle='', linewidth=0, marker='o', fillstyle='none', mew=0.5)

        ax.plot(profile['x']/h,profile['z']/h,color="#222222")
        ax.fill_between(profile['x']/h,profile['z']/h, color="#222222", alpha=0.8)
        
        # 1. Extract all handles and labels from the current axis
        handles, labels = ax.get_legend_handles_labels()

        # 2. Use a dict comprehension to retain only unique labels
        unique_legend = dict(zip(labels, handles))

        # 3. Pass unique elements back into the legend
        ax.legend(unique_legend.values(), unique_legend.keys())

        plt.axis('equal')
        
        if pt == "h_mean":
            plt.title("Scaled Time and Horizontally-Averaged Velocity")
        if pt == "center_slice":
            plt.title("Scaled Time-Averaged Velocity")

        plt.xlabel("x/h with Mean Velocity (U/U_b)")
        plt.ylabel("z/h")
        plt.savefig('ufr_' + pt + '_flow.png')


def calc_error(output, expdata):
    """
    Compare the model spanwise-mean velocity profiles in ``output[name]['h_mean']``
    to the experimental profiles in ``expdata``.

    A nearest-offset match is used to pair each experimental x/h station with the
    closest model offset in physical space.
    """
    h = 0.028
    errors = {}

    for name, m in output.items():
        case_errors = []
        u_b = m['ummid']/1.0593

        for o in range(m['noffsets']):
            model_x = m['offsets'][o] / h
            model_profile = np.asarray(m['h_mean'][o], dtype=float)
            model_z = np.asarray(m['zvals'], dtype=float) / h


            # Find the closest experimental station in x/h
            exp_x = min(expdata.keys(), key=lambda key: abs(float(key) - model_x*h))
            exp_df = expdata[exp_x]

            exp_z = np.asarray(exp_df[0], dtype=float)
            exp_u = np.asarray(exp_df[1], dtype=float)

            # Interpolate the model profile onto the experimental z/h grid
            model_interp = np.interp(exp_z, model_z, model_profile)


            model_scaled = model_interp / (u_b)
            exp_scaled = exp_u 

            rmse = np.sqrt(np.mean((model_scaled - exp_scaled) ** 2))
            mae = np.mean(np.abs(model_scaled - exp_scaled))
            max_abs = np.max(np.abs(model_scaled - exp_scaled))

            case_errors.append({
                'offset': m['offsets'][o],
                'x_h': model_x,
                'exp_x_h': float(exp_x/h),
                'rmse': rmse,
                'mae': mae,
                'max_abs': max_abs,
                'case_label' : m['label']
            })

        errors[name] = case_errors

        if case_errors:
            mean_rmse = np.mean([e['rmse'] for e in case_errors])
            mean_mae = np.mean([e['mae'] for e in case_errors])
            print(f"{name}: mean RMSE={mean_rmse:.6g}, mean MAE={mean_mae:.6g}")

    return errors


def bar_plot_errors(errors):
    """
    Create one bar plot per error metric.

    Each plot uses the offset location on the x-axis, and every offset contains a
    cluster of bars for the available case models. The bar color identifies the
    case/model, while the per-plot figure corresponds to a single error type.
    """
    if not errors:
        return

    metric_names = ['rmse', 'mae', 'max_abs']
    case_names = list(errors.keys())
    offsets = sorted({float(entry['x_h']) for case_errors in errors.values() for entry in case_errors})

    case_colors = {
        case_name: color
        for case_name, color in zip(
            case_names,
            [
                '#D55E00',
                '#56B4E9',
                '#CC79A7',
                '#009E73',
                '#F0E442',
                '#0072B2',
                '#E69F00',
                '#999999',
            ],
        )
    }

    for metric_name in metric_names:
        fig, ax = plt.subplots(figsize=(12, 8))
        x = np.arange(len(offsets))
        ncases = len(case_names)
        cluster_width = 0.72 / max(1, ncases)

        for case_index, case_name in enumerate(case_names):
            color = case_colors.get(case_name, '#777777')
            case_map = {float(entry['x_h']): entry for entry in errors[case_name]}
            values = [case_map.get(offset, {}).get(metric_name, np.nan) for offset in offsets]
            positions = x + (case_index - (ncases - 1) / 2.0) * cluster_width

            ax.bar(
                positions,
                values,
                width=cluster_width,
                color=color,
                edgecolor='black',
                linewidth=0.5,
                alpha=0.75,
                label=errors[case_name][0]['case_label'],
            )

        ax.set_xticks(x)
        ax.set_xticklabels([f"{offset:.2f}" for offset in offsets])
        ax.set_xlabel('Offset location (x/h)')
        ax.set_ylabel(metric_name.upper() + ' (m/s)')
        ax.set_title(f'{metric_name.upper()} by offset location')
        ax.grid(axis='y', linestyle='--', alpha=0.4)
        ax.legend(loc='upper right', fontsize=10)

        plt.tight_layout()
        plt.savefig(f'error_bar_plot_{metric_name}.png')
        plt.close(fig)


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

            if c['type'] == 'sgf':
                output[c['name']] = make_means_sgf(planedata[c['name']],c['name'],c['opt'],c['label'])
            else:
                output[c['name']] = make_means_ugf(planedata[c['name']],c['name'],c['opt'],c['label'])

            print("Plane output steps:",output[c['name']]['ns'])

            if c['type'] == 'ugf':
                force_data = pd.read_csv(c['path'] + '/forcing.dat', sep=r'\s+', skipinitialspace=True)
                print("Mean Forcing:",force_data.iloc[-1000:]['Fx'].mean()/1000.0)

    plot_planes_all(output,expdata,profile)
    errors = calc_error(output,expdata)
    bar_plot_errors(errors)
        

if __name__=="__main__":
    main()