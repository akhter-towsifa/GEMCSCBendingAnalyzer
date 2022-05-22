import os

#filepath = "/eos/user/d/daebi/cosmic_MC/" #MC cosmic
#filepath = "/eos/cms/store/express/Commissioning2020/ExpressCosmics/FEVT/Express-v1/000/337/973/00000/" #MWGR4 express dataset
#filepath = "/eos/cms/store/express/Commissioning2020/ExpressCosmics/FEVT/Express-v1/000/338/714/" #MWGR5 express dataset
#filepath = "/eos/cms/store/express/Commissioning2021/ExpressCosmics/FEVT/Express-v1/000/341/343/" #MWGR3 express dataset
filepath = "/eos/cms/store/express/Commissioning2022/ExpressCosmics/FEVT/Express-v1/000/" #CRAFT
runlist = [348773,348776,348777,348838,348908,348955,349016,349073,349078,349079,349084,349146,349147,349260,349263,349348,349422,349433,349435,349436,349437,349527,349528,349611,349703,349758,349833,349834,349839,349840,349893,349963,350010,350107,350142,350166,350174,350254,350294,350361,350424,350425,350431,350450,350459,350460,350462,350463,350464,350490,350491,350561,350619]
pwd = os.getcwd()+'/'

submit_together = True #cd scripts // condor_submit submit_all_condor.sh
submit_individ = False #./submit_all.sh
n_files_per_script = 15
jobname = "ME11Iter2_GE11Iter0_210522"
high_priority = True

#Whether or not to use a misalignment
misalign = False
misalign_db = "gemAl.db"
misalign_tag = "GEM"
misalign_APR_tag = "test"

#Set up the current folder to be the working space
curr_dir = pwd + jobname + "/"
os.system("mkdir -p {curr_dir}".format(curr_dir = curr_dir))

sub_all_name = "{jobname}/submit_all.sh".format(jobname = jobname)
sub_all = open(sub_all_name, "w")
sub_all.write("#!/bin/bash\n")
sub_all.write("date\n")
sub_all.write("cd "+curr_dir+"\n")
sub_all.write("eval `scramv1 runtime -sh`\n")

ntotal_jobs = 0

for run in runlist:
  job_filelist = []
  files_in_script = 0
  runjob_counter = 0

  run_str = str(run)
  run_inputdir = filepath + run_str[:3] + "/" + run_str[3:] + "/"
  run_dir = curr_dir + run_str + "/"
  os.system("mkdir -p {run_dir}".format(run_dir = run_dir))
  os.system("mkdir -p {run_dir}/output".format(run_dir = run_dir))
  os.system("mkdir -p {run_dir}/error".format(run_dir = run_dir))
  os.system("mkdir -p {run_dir}/log".format(run_dir = run_dir))
  os.system("mkdir -p {run_dir}/scripts".format(run_dir = run_dir))
  os.system("mkdir -p {run_dir}/rootfiles".format(run_dir = run_dir))

  for subdir in os.listdir(run_inputdir):
    sub_run_inputdir = run_inputdir + subdir + "/"
    for tmp_filename in os.listdir(sub_run_inputdir):
      if ".root" not in tmp_filename: continue
      filename = sub_run_inputdir + tmp_filename

      if files_in_script == 0:
        scriptname = run_dir+"scripts/run"+run_str+"_job{njob}.sh".format(njob = runjob_counter)
        c_script = open(scriptname, "w")
        runjob_counter += 1

      job_filelist.append('file:'+filename)
      files_in_script += 1

      if files_in_script == n_files_per_script or (subdir == os.listdir(run_inputdir)[-1] and tmp_filename == os.listdir(sub_run_inputdir)[-1]):
        #print("Got filelist, writing condor script for {nfiles} files".format(nfiles = files_in_script))
        files_in_script = 0

        outname = run_dir+"rootfiles/out_run"+run_str+"_job{njob}.root".format(njob = runjob_counter)
        filelist = ""
        for fname in job_filelist:
          filelist += "file:" + fname
          if fname != job_filelist[-1]:
            filelist += ","
        c_script.write("#!/bin/bash\n")
        c_script.write("date\n")
        c_script.write("cd "+pwd+"\n")
        c_script.write("eval `scramv1 runtime -sh`\n")
        c_script.write("cmsRun run_ME11_GE11_condor.py inputFiles={filelist} outputFile={outname}\n".format(filelist = filelist, outname = outname))
        c_script.write("echo 'done'\n")
        c_script.write("date\n")

        job_filelist = []
        os.system("chmod 755 "+scriptname)

        ntotal_jobs += 1
  print("Run " + run_str + " had {njob} total jobs".format(njob = runjob_counter))

  submit_run_name = "{run_dir}/scripts/submit_run{run_str}.sh".format(run_dir = run_dir, run_str = run_str)
  submit_run = open(submit_run_name, "w")
  if high_priority:
    submit_run.write("""executable              = $(filename)
output                  = {run_dir}/output/$(filename).out
error                   = {run_dir}/error/$(filename).err
log                     = {run_dir}/log/$(filename).log
request_memory          = 4000M
+JobFlavour             = "workday"
universe                = vanilla
+AccountingGroup        = "group_u_CMS.CAF.ALCA"
queue filename matching files run*.sh""".format(run_dir = run_dir))

  else:
    submit_run.write("""executable              = $(filename)
output                  = {run_dir}/output/$(filename).out
error                   = {run_dir}/error/$(filename).err
log                     = {run_dir}/log/$(filename).log
request_memory          = 4000M
+JobFlavour             = "workday"
universe                = vanilla
queue filename matching files run*.sh""".format(run_dir = run_dir))

  os.system("chmod 755 "+submit_run_name)

  sub_all.write("echo run {run_str}\n".format(run_str = run_str))
  sub_all.write("cd {run_dir}/scripts/\n".format(run_dir = run_dir))
  sub_all.write("condor_submit {run_dir}/scripts/submit_run{run_str}.sh\n".format(run_dir = run_dir, run_str = run_str))

print("Total jobs made was ", ntotal_jobs)

sub_all.write("echo 'done'\n")
sub_all.write("cd "+curr_dir+"\n")
sub_all.write("date\n")

os.system("chmod 755 "+sub_all_name)

print("Total number of runs was ", len(runlist))
