
pattern = '...........((((((((((((((((....&)))))).)))))).)))).'

target = pattern.split("&")[0]
smRNA = pattern.split("&")[1]
smRNA_rev = smRNA[::-1]

smRNA_len = len(smRNA_rev)
sub_target = target[-smRNA_len:]

print(sub_target)
print(smRNA)

