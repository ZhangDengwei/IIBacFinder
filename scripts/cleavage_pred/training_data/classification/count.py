import json

with open('/data2/ygao/metatrain3/data.json', 'r') as f:
    data = json.load(f)
    _class={}
    for dic in data:
        if dic['label'] in _class:
           _class['%s'%dic['label']]+=1
        else:
           _class['%s'%dic['label']]=1
    print(_class)
