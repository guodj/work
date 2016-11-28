import gitm
import gitm_3D_global_plots as gpt

a = gitm.GitmBin('/home/gdj/tmp/3DALL_t100516_113000.bin')
gpt.gitm_single_3D_image('polar', 'Temperature', a)
