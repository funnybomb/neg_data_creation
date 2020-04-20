#!/usr/bin/env python
# coding: utf-8

# In[1]:


from Bio import pairwise2 as pw
from Bio import Seq
from Bio import SeqIO 
import Bio
from Bio.Alphabet import IUPAC
import numpy as np
import pandas as pd
import re


# In[2]:


import time
from functools import wraps
  
def fn_timer(function):
  @wraps(function)
  def function_timer(*args, **kwargs):
    t0 = time.time()
    result = function(*args, **kwargs)
    t1 = time.time()
    print ("Total time running %s: %s seconds" %
        (function.__name__, str(t1-t0))
        )
    return result
  return function_timer


# In[3]:


uniprot = list(SeqIO.parse('uniprot-proteome_human.fasta','fasta'))


# In[4]:


for seq in uniprot:   #
    x = seq.id
    seq.id = x.split("|")[1]


# In[5]:


ids = [seq.id for seq in uniprot]
names = [seq.name for seq in uniprot]


# In[6]:


uni_dict = {}
for i in uniprot:
    uni_dict[i.id] =i


# In[7]:


grist = pd.read_table('gristone_positive_data.txt')


# In[8]:


grist.head(3)


# In[9]:


#@fn_timer
def calc_identity_score(seq1,seq2,gap=-0.5,extend=-0.1):
    """
    return seq1 seq1 pairwise identiy score
    seq1 is positive 
    seq1 seq2 is string

    Bio.pairwise2.format_alignment output:
    MPKGKKAKG------
      |||||||      
    --KGKKAKGKKVAPA
      Score=7
      
    alignment output: [('MPKGKKAKG------', '--KGKKAKGKKVAPA', 7.0, 0, 15)]
    score = ali[0][2] = 7
    """
    ali = pw.align.globalxs(seq1,seq2,gap,extend,score_only=True)
    # gap penalty = -0.5 in case of cak caak score =3
    return ali/min(len(seq1),len(seq2)) # 返回短序列的值,防止substring


# In[10]:


#@fn_timer
def window_generator(seq, window_lenth=7,step=1):
    """
    return list of seq window slide
    seq is string
    """
    if len(seq) >= window_lenth:
        return [seq[i:i+window_lenth] for i in range(0,len(seq)-window_lenth+1,step)]
    else:
        return []


# In[11]:


window_generator('123',3)


# In[12]:


from Bio.pairwise2 import format_alignment


# In[13]:


def flat(nums):
    res = []
    for i in nums:
        if isinstance(i, list):
            res.extend(flat(i))
        else:
            res.append(i)
    return res


# In[14]:


#@fn_timer
def slide_with_flank(seq,full,step=1,up_flank=6,down_flank=6,flags=0):
    """
    return window slide result as list for a full str given potential sub seq and it's flank removed
    seq=abc, full = 01234abc45678abc1234 up=1 down=1 step =1
    result is ['012', '123', '567', '234']
    生成一个str full 的滑动窗口自片段,去除了seq和其上下游的部分
    """
    res = []
    window_len = len(seq)
    coords = [i.span() for i in re.finditer(seq,full,flags)] # = search_all
    # 处理首尾情况
    if len(coords) == 0:
        res.append(window_generator(full,window_lenth=window_len,step=step))

    elif len(coords) == 1:
        if (coords[0][0]-up_flank) >= 0:
            res.append(window_generator(full[0:coords[0][0]-up_flank],
                                        window_lenth=window_len,step=step))
        if (coords[0][1]+down_flank) <= len(full):
            res.append(window_generator(full[coords[0][1]+up_flank:],
                                        window_lenth=window_len,step=step))
    else: # len(coords) >1  
        if (coords[0][0]-up_flank) >= 0:
            res.append(window_generator(full[0:coords[0][0]-up_flank],
                                        window_lenth=window_len,step=step))
        for i in range(1,len(coords)):
        ## 处理 1 2 之间的东西
            if coords[i][0] - coords[i-1][1] > up_flank+down_flank+window_len:
                res.append(window_generator(full[coords[i-1][1]+up_flank:coords[i][0]-down_flank],
                                            window_lenth=window_len,step=step))
                           
        if (coords[-1][1]+down_flank) <= len(full):
            res.append(window_generator(full[coords[-1][1]+down_flank:],
                                        window_lenth=window_len,step=step))           
        
    return flat(res)


# In[15]:


slide_with_flank('abc','01234abc456raa78abc7779',up_flank=2,down_flank=2)


# In[16]:


###### @fn_timer
def filter_with_identity_affinity(seq,full,identity_cutoff=0.5,step=1,up=6,down=6,flags=0):
    filtered = {}
    slides = slide_with_flank(seq,full,step=step,up_flank=up,down_flank=down,flags=flags)
    for s in slides:
        identity_score = calc_identity_score(seq,s)
        if identity_score <= identity_cutoff:
            #if s in filtered.keys(): #已经有score, 寸大值表明有大片段
            filtered[s] = identity_score
    return filtered


# In[ ]:


grist['neg1'] =np.full_like(grist['peptide'],'o')
t0 = time.process_time()
for i in range(0,len(grist)):
#for i in range(0,10000):
    g = grist.iloc[i]
    if g['uniport_id'] in ids:
        pep = str(g['left_flanking']) + str(g['peptide']) + str(g['right_flanking'])
        full = str(uni_dict[g['uniport_id']].seq)
        res_local = filter_with_identity_affinity(pep,full,step=10)
#        res_full = {}
#        for win in res_local.keys(): # 过滤掉阳性集在window结果中 identity>0.5的
#            for pep in np.random.choice(grist['peptide'],100):#.remove(g['peptide']):
#                score = calc_identity_score(win,pep)
#                if score <0.5:
#                    res_full[win] = score
#                else:
#                    break
#        res_full_key = sorted(res_full.keys(),key=lambda x:(res_local[x],res_full[x]))
        k=sorted(res_local.keys(),key=lambda x:res_local[x])
        if len(k):
            grist.loc[i,'neg1'] = k[0] 
        
print(time.process_time()-t0)
#grist[grist['neg1'] != 'o'].to_csv('testdata0401.txt',index = False)


# In[18]:


grist[['left_flanking','right_flanking']]


# In[28]:


a=grist['neg1'][0]
#gri


# In[30]:


get_ipython().run_line_magic('pinfo', 'grist.loc')


# In[25]:


a


# In[ ]:





# In[22]:


get_ipython().run_line_magic('pinfo', 'grist.loc')


# In[53]:


grist.head()


# In[22]:



persons={'ZhangSan':'male',
         'LiSi':'male',
         'WangHong':'female'}

#找出所有男性
males = filter(lambda x:'male'== x[1], persons.items())

for (key,value) in males:
    print('%s : %s' % (key,value))


# In[42]:


get_ipython().run_line_magic('pinfo', 'grist.to_csv')


# In[2]:


get_ipython().run_cell_magic('writefile', 'train.py', '')


# In[ ]:




