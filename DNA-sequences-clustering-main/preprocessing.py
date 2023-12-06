with open("DNA-sequences-clustering-main\DNA-sequences-clustering-main\data files\human_short_unlabeled.txt", "r") as f:
    lines = f.readlines()
    trimmed_lines = [line.split('\t', 1)[0] + '\n' for line in lines]
   
with open("DNA-sequences-clustering-main\DNA-sequences-clustering-main\data files\human_short_unlabeled.txt", "w") as f:
    f.writelines(trimmed_lines)