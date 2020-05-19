import os

if __name__ == '__main__':

    # get the complete directory list
    complete_dir_list = os.listdir()

    # get the dat file list
    dat_file_list = []
    for i, val in enumerate(complete_dir_list):
        if ".dat" in val:
            dat_file_list.append(val)

    min_sys_size = 100000000
    max_sys_size = 0
    for i, val in enumerate(dat_file_list):
        with open(val, 'r') as f:
            lines = f.read().splitlines()
            last_line = lines[-1]
            if "Total" not in last_line:
                print(val)
                # print(last_line)
                debased_val = str(val.replace("log_ground_state_", "").split(".dat", 1)[0])
                debased_val_entries = debased_val.split('_')
                q = int(debased_val_entries[debased_val_entries.index("nphi") + 2])
                LxMUC = int(debased_val_entries[debased_val_entries.index("LxMUC") + 1])
                Ly = int(debased_val_entries[debased_val_entries.index("Ly") + 1])
                sys_size = q*LxMUC*Ly
                print("System size =", sys_size)
                if sys_size < min_sys_size:
                    min_sys_size = sys_size
                elif sys_size > max_sys_size:
                    max_sys_size = sys_size
    print("Min/max system size = (", min_sys_size, ",", max_sys_size, ")")
