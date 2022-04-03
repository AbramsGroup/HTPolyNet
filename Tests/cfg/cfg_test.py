import HTPolyNet.configuration as c

print('JSON FORMAT----------------')
C=c.Configuration.read('test.json')
print(C)
print('YAML FORMAT----------------')
C=c.Configuration.read('test.yaml')
print(C)
print('TXT FORMAT----------------')
C=c.Configuration.read('test.txt')
print(C)