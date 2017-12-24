
Pymaceutical exercise data analysis:
-By analysing the Tumor volume (mm3) variation with the different drugs during the 45 days, we can verify that the drug Capomulin caused a decrease of about 20% in the Tumor volume
-The spread of  Metastatic sites was also lower with the Capomulin drug over the 45 days
-The survival rate with treatment with drug Capomulin was roughly of 80%, whereas it ranged between 40-60% with the other treatments

Overall we can conclude that Capomulin was the drug that had better effects in terms of reducing the size of the volume of the tumors, reducing the spreading of metastatic sites and increasing the values of survival rates of mice.


```python
#Import Dependencies
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os
```


```python
# Take in all of raw_data and read it into pandas
clinical_trial=os.path.join('raw_data', 'clinicaltrial_data.csv')
mouse_drug_data=os.path.join('raw_data', 'mouse_drug_data.csv')
clinical_trial_df = pd.read_csv(clinical_trial)
mouse_drug_data_df = pd.read_csv(mouse_drug_data)
clinical_trial_df.head()
mouse_drug_data_df.head()
```




<div>
<style>
    .dataframe thead tr:only-child th {
        text-align: right;
    }

    .dataframe thead th {
        text-align: left;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>Mouse ID</th>
      <th>Drug</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>f234</td>
      <td>Stelasyn</td>
    </tr>
    <tr>
      <th>1</th>
      <td>x402</td>
      <td>Stelasyn</td>
    </tr>
    <tr>
      <th>2</th>
      <td>a492</td>
      <td>Stelasyn</td>
    </tr>
    <tr>
      <th>3</th>
      <td>w540</td>
      <td>Stelasyn</td>
    </tr>
    <tr>
      <th>4</th>
      <td>v764</td>
      <td>Stelasyn</td>
    </tr>
  </tbody>
</table>
</div>




```python
#Merge data of both files and format
clinical_mouse_merge_df = pd.merge(clinical_trial_df, mouse_drug_data_df, on="Mouse ID", how="outer")
clinical_mouse_merge_df['Tumor Volume (mm3)'] = clinical_mouse_merge_df['Tumor Volume (mm3)'].map("{:.1f}".format)
#clinical_mouse_merge_df["Drug"].unique()
clinical_mouse_merge_df["Tumor Volume (mm3)"]= pd.to_numeric(clinical_mouse_merge_df["Tumor Volume (mm3)"])
clinical_mouse_merge_df["Timepoint"]= pd.to_numeric(clinical_mouse_merge_df["Timepoint"])
clinical_mouse_merge_df["Metastatic Sites"]= pd.to_numeric(clinical_mouse_merge_df["Metastatic Sites"])
clinical_mouse_merge_df.sort_values("Timepoint").head()
```




<div>
<style>
    .dataframe thead tr:only-child th {
        text-align: right;
    }

    .dataframe thead th {
        text-align: left;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>Mouse ID</th>
      <th>Timepoint</th>
      <th>Tumor Volume (mm3)</th>
      <th>Metastatic Sites</th>
      <th>Drug</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>b128</td>
      <td>0</td>
      <td>45.0</td>
      <td>0</td>
      <td>Capomulin</td>
    </tr>
    <tr>
      <th>1535</th>
      <td>i635</td>
      <td>0</td>
      <td>45.0</td>
      <td>0</td>
      <td>Propriva</td>
    </tr>
    <tr>
      <th>565</th>
      <td>g791</td>
      <td>0</td>
      <td>45.0</td>
      <td>0</td>
      <td>Ramicane</td>
    </tr>
    <tr>
      <th>1545</th>
      <td>w746</td>
      <td>0</td>
      <td>45.0</td>
      <td>0</td>
      <td>Propriva</td>
    </tr>
    <tr>
      <th>1547</th>
      <td>r107</td>
      <td>0</td>
      <td>45.0</td>
      <td>0</td>
      <td>Propriva</td>
    </tr>
  </tbody>
</table>
</div>




```python
#Tumor Response to Treatment, mean of the Tumor volume (mm3)
Drug_effect_mean= clinical_mouse_merge_df.groupby(["Drug","Timepoint"]).mean()
Drug_effect_mean_df=pd.DataFrame(Drug_effect_mean['Tumor Volume (mm3)'])
Drug_effect_mean_df.head()
```




<div>
<style>
    .dataframe thead tr:only-child th {
        text-align: right;
    }

    .dataframe thead th {
        text-align: left;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th></th>
      <th>Tumor Volume (mm3)</th>
    </tr>
    <tr>
      <th>Drug</th>
      <th>Timepoint</th>
      <th></th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th rowspan="5" valign="top">Capomulin</th>
      <th>0</th>
      <td>45.000000</td>
    </tr>
    <tr>
      <th>5</th>
      <td>44.268000</td>
    </tr>
    <tr>
      <th>10</th>
      <td>43.080000</td>
    </tr>
    <tr>
      <th>15</th>
      <td>42.066667</td>
    </tr>
    <tr>
      <th>20</th>
      <td>40.717391</td>
    </tr>
  </tbody>
</table>
</div>




```python
#Tumor Response to Treatment, stardard error of the mean of the Tumor volume (mm3)
Drug_effect_sem= clinical_mouse_merge_df.groupby(["Drug","Timepoint"]).sem()
Drug_effect_sem_df=pd.DataFrame(Drug_effect_sem['Tumor Volume (mm3)'])
Drug_effect_sem_df.head()
```




<div>
<style>
    .dataframe thead tr:only-child th {
        text-align: right;
    }

    .dataframe thead th {
        text-align: left;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th></th>
      <th>Tumor Volume (mm3)</th>
    </tr>
    <tr>
      <th>Drug</th>
      <th>Timepoint</th>
      <th></th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th rowspan="5" valign="top">Capomulin</th>
      <th>0</th>
      <td>0.000000</td>
    </tr>
    <tr>
      <th>5</th>
      <td>0.449137</td>
    </tr>
    <tr>
      <th>10</th>
      <td>0.704036</td>
    </tr>
    <tr>
      <th>15</th>
      <td>0.837107</td>
    </tr>
    <tr>
      <th>20</th>
      <td>0.908820</td>
    </tr>
  </tbody>
</table>
</div>




```python
#Create pivor table of the Tumor Response to Treatment, mean of the Tumor volume (mm3)
Drug_effect_tumor_mean_pivot_table=Drug_effect_mean_df.pivot_table('Tumor Volume (mm3)', ['Timepoint'], 'Drug')
Drug_effect_tumor_mean_pivot_table.columns
```




    Index(['Capomulin', 'Ceftamin', 'Infubinol', 'Ketapril', 'Naftisol', 'Placebo',
           'Propriva', 'Ramicane', 'Stelasyn', 'Zoniferol'],
          dtype='object', name='Drug')




```python
#Create pivor table of the Tumor Response to Treatment, stardard error of the mean of the Tumor volume (mm3) 
Drug_effect_tumor_sem_pivot_table=Drug_effect_sem_df.pivot_table('Tumor Volume (mm3)', ['Timepoint'], 'Drug')
Drug_effect_tumor_sem_pivot_table.head()
```




<div>
<style>
    .dataframe thead tr:only-child th {
        text-align: right;
    }

    .dataframe thead th {
        text-align: left;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>Drug</th>
      <th>Capomulin</th>
      <th>Ceftamin</th>
      <th>Infubinol</th>
      <th>Ketapril</th>
      <th>Naftisol</th>
      <th>Placebo</th>
      <th>Propriva</th>
      <th>Ramicane</th>
      <th>Stelasyn</th>
      <th>Zoniferol</th>
    </tr>
    <tr>
      <th>Timepoint</th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
    </tr>
    <tr>
      <th>5</th>
      <td>0.449137</td>
      <td>0.164096</td>
      <td>0.235533</td>
      <td>0.266584</td>
      <td>0.201272</td>
      <td>0.216834</td>
      <td>0.233172</td>
      <td>0.485356</td>
      <td>0.240175</td>
      <td>0.189584</td>
    </tr>
    <tr>
      <th>10</th>
      <td>0.704036</td>
      <td>0.234734</td>
      <td>0.283339</td>
      <td>0.356324</td>
      <td>0.321307</td>
      <td>0.400980</td>
      <td>0.377399</td>
      <td>0.719240</td>
      <td>0.433280</td>
      <td>0.263944</td>
    </tr>
    <tr>
      <th>15</th>
      <td>0.837107</td>
      <td>0.334288</td>
      <td>0.352805</td>
      <td>0.577249</td>
      <td>0.444669</td>
      <td>0.612733</td>
      <td>0.466272</td>
      <td>0.770151</td>
      <td>0.493940</td>
      <td>0.367334</td>
    </tr>
    <tr>
      <th>20</th>
      <td>0.908820</td>
      <td>0.356109</td>
      <td>0.476416</td>
      <td>0.725309</td>
      <td>0.594510</td>
      <td>0.839567</td>
      <td>0.554597</td>
      <td>0.787014</td>
      <td>0.619751</td>
      <td>0.533493</td>
    </tr>
  </tbody>
</table>
</div>




```python
#Comparison analysis for ["Capomulin","Infubinol","Ketapril","Placebo" ]
#mean
Comp_mean_Drug_effect=Drug_effect_tumor_mean_pivot_table[["Capomulin","Infubinol","Ketapril","Placebo" ]]
#sem
Comp_sem_Drug_effect=Drug_effect_tumor_sem_pivot_table[["Capomulin","Infubinol","Ketapril","Placebo" ]]
Comp_sem_Drug_effect ['Capomulin']
```




    Timepoint
    0     0.000000
    5     0.449137
    10    0.704036
    15    0.837107
    20    0.908820
    25    0.884443
    30    0.933820
    35    1.052804
    40    1.223684
    45    1.226053
    Name: Capomulin, dtype: float64




```python
#def a function to make a scatter plot for one of the drugs

def scatterplot(x_data, y_data, error_data, color, marker):
    plt.errorbar(x_data, y_data,yerr = error_data,capsize=5, capthick=2,ecolor='k'
                 , color=colors[color_index],marker=markers[marker_index]
                 ,markersize=10,linestyle='--', dashes=(3, 1),linewidth=2,alpha = 0.75)
                
colors = ['g','r','b','y']
markers=["^","s","d","o"]
color_index = 0
marker_index=0
#xlabel="Time (Days)"
scatterplot(x_data = Comp_mean_Drug_effect.index
            , y_data = Comp_mean_Drug_effect['Capomulin']
          , error_data = Comp_sem_Drug_effect ['Capomulin']
            ,color = colors[color_index], marker=markers[marker_index])  
color_index += 1
marker_index +=1
plt.title("Tumor Response to Treatment")
plt.xlabel("Time (Days)")
plt.ylabel("Tumor Volume (mm3)")
plt.legend(loc="best")
plt.grid(True, linestyle='dotted')
plt.xticks(range(0,50, 5))
plt.yticks(range(20,90, 10))
plt.xlim(xmin=0)

plt.show()
```


![png](output_9_0.png)



```python
#Make a scatter plot with comparison analysis for ["Capomulin","Infubinol","Ketapril","Placebo" ] 
#based on the function above
colors = ['g','r','b','y']
markers=["^","s","d","o"]
color_index = 0
marker_index=0
for drug in Comp_mean_Drug_effect:
    scatterplot(x_data = Comp_mean_Drug_effect.index, 
                y_data = Comp_mean_Drug_effect[drug], error_data = Comp_sem_Drug_effect[drug],color = colors[color_index],marker=markers[marker_index])
    color_index += 1
    marker_index +=1
    
xmin= Comp_mean_Drug_effect.index.min()
xmax=Comp_mean_Drug_effect.index.max()
plt.title("Tumor Response to Treatment")
plt.xlabel("Time (Days)")
plt.ylabel("Tumor Volume (mm3)")
plt.legend(loc="best")
plt.grid(True, linestyle='dotted')
plt.xticks(range(0,50, 5))
plt.yticks(range(20,90, 10))
plt.xlim(xmin, xmax)
plt.savefig("Fig1.png")

plt.show() 

```


![png](output_10_0.png)



```python
#Metastatic Sites Response to Treatment, mean of the Metastatic Sites
Drug_effect_mean= clinical_mouse_merge_df.groupby(["Drug","Timepoint"]).mean()
Metastatic_effect_mean_df=pd.DataFrame(Drug_effect_mean['Metastatic Sites'])
Metastatic_effect_mean_df.head()

```




<div>
<style>
    .dataframe thead tr:only-child th {
        text-align: right;
    }

    .dataframe thead th {
        text-align: left;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th></th>
      <th>Metastatic Sites</th>
    </tr>
    <tr>
      <th>Drug</th>
      <th>Timepoint</th>
      <th></th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th rowspan="5" valign="top">Capomulin</th>
      <th>0</th>
      <td>0.000000</td>
    </tr>
    <tr>
      <th>5</th>
      <td>0.160000</td>
    </tr>
    <tr>
      <th>10</th>
      <td>0.320000</td>
    </tr>
    <tr>
      <th>15</th>
      <td>0.375000</td>
    </tr>
    <tr>
      <th>20</th>
      <td>0.652174</td>
    </tr>
  </tbody>
</table>
</div>




```python
#Metastatic Sites Response, stardard error of the mean of Metastatic Sites
Drug_effect_sem= clinical_mouse_merge_df.groupby(["Drug","Timepoint"]).sem()
Metastatic_effect_sem_df=pd.DataFrame(Drug_effect_sem['Metastatic Sites'])
Metastatic_effect_sem_df.head()
```




<div>
<style>
    .dataframe thead tr:only-child th {
        text-align: right;
    }

    .dataframe thead th {
        text-align: left;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th></th>
      <th>Metastatic Sites</th>
    </tr>
    <tr>
      <th>Drug</th>
      <th>Timepoint</th>
      <th></th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th rowspan="5" valign="top">Capomulin</th>
      <th>0</th>
      <td>0.000000</td>
    </tr>
    <tr>
      <th>5</th>
      <td>0.074833</td>
    </tr>
    <tr>
      <th>10</th>
      <td>0.125433</td>
    </tr>
    <tr>
      <th>15</th>
      <td>0.132048</td>
    </tr>
    <tr>
      <th>20</th>
      <td>0.161621</td>
    </tr>
  </tbody>
</table>
</div>




```python
#Create pivor table of the Metastatic Sites Response to Treatment, mean of the Metastatic Sites
Metastatic_effect_mean_pivot_table=Metastatic_effect_mean_df.pivot_table('Metastatic Sites', ['Timepoint'], 'Drug')
Metastatic_effect_mean_pivot_table.head()
```




<div>
<style>
    .dataframe thead tr:only-child th {
        text-align: right;
    }

    .dataframe thead th {
        text-align: left;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>Drug</th>
      <th>Capomulin</th>
      <th>Ceftamin</th>
      <th>Infubinol</th>
      <th>Ketapril</th>
      <th>Naftisol</th>
      <th>Placebo</th>
      <th>Propriva</th>
      <th>Ramicane</th>
      <th>Stelasyn</th>
      <th>Zoniferol</th>
    </tr>
    <tr>
      <th>Timepoint</th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
    </tr>
    <tr>
      <th>5</th>
      <td>0.160000</td>
      <td>0.380952</td>
      <td>0.280000</td>
      <td>0.304348</td>
      <td>0.260870</td>
      <td>0.375000</td>
      <td>0.320000</td>
      <td>0.120000</td>
      <td>0.240000</td>
      <td>0.166667</td>
    </tr>
    <tr>
      <th>10</th>
      <td>0.320000</td>
      <td>0.600000</td>
      <td>0.666667</td>
      <td>0.590909</td>
      <td>0.523810</td>
      <td>0.833333</td>
      <td>0.565217</td>
      <td>0.250000</td>
      <td>0.478261</td>
      <td>0.500000</td>
    </tr>
    <tr>
      <th>15</th>
      <td>0.375000</td>
      <td>0.789474</td>
      <td>0.904762</td>
      <td>0.842105</td>
      <td>0.857143</td>
      <td>1.250000</td>
      <td>0.764706</td>
      <td>0.333333</td>
      <td>0.782609</td>
      <td>0.809524</td>
    </tr>
    <tr>
      <th>20</th>
      <td>0.652174</td>
      <td>1.111111</td>
      <td>1.050000</td>
      <td>1.210526</td>
      <td>1.150000</td>
      <td>1.526316</td>
      <td>1.000000</td>
      <td>0.347826</td>
      <td>0.952381</td>
      <td>1.294118</td>
    </tr>
  </tbody>
</table>
</div>




```python
#Create pivor table of the  Metastatic Sites to Response Treatment, stardard error of the mean of the Metastatic Sites
Metastatic_effect_sem_pivot_table=Metastatic_effect_sem_df.pivot_table('Metastatic Sites', ['Timepoint'], 'Drug')
Metastatic_effect_sem_pivot_table.head()
```




<div>
<style>
    .dataframe thead tr:only-child th {
        text-align: right;
    }

    .dataframe thead th {
        text-align: left;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>Drug</th>
      <th>Capomulin</th>
      <th>Ceftamin</th>
      <th>Infubinol</th>
      <th>Ketapril</th>
      <th>Naftisol</th>
      <th>Placebo</th>
      <th>Propriva</th>
      <th>Ramicane</th>
      <th>Stelasyn</th>
      <th>Zoniferol</th>
    </tr>
    <tr>
      <th>Timepoint</th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
    </tr>
    <tr>
      <th>5</th>
      <td>0.074833</td>
      <td>0.108588</td>
      <td>0.091652</td>
      <td>0.098100</td>
      <td>0.093618</td>
      <td>0.100947</td>
      <td>0.095219</td>
      <td>0.066332</td>
      <td>0.087178</td>
      <td>0.077709</td>
    </tr>
    <tr>
      <th>10</th>
      <td>0.125433</td>
      <td>0.152177</td>
      <td>0.159364</td>
      <td>0.142018</td>
      <td>0.163577</td>
      <td>0.115261</td>
      <td>0.105690</td>
      <td>0.090289</td>
      <td>0.123672</td>
      <td>0.109109</td>
    </tr>
    <tr>
      <th>15</th>
      <td>0.132048</td>
      <td>0.180625</td>
      <td>0.194015</td>
      <td>0.191381</td>
      <td>0.158651</td>
      <td>0.190221</td>
      <td>0.136377</td>
      <td>0.115261</td>
      <td>0.153439</td>
      <td>0.111677</td>
    </tr>
    <tr>
      <th>20</th>
      <td>0.161621</td>
      <td>0.241034</td>
      <td>0.234801</td>
      <td>0.236680</td>
      <td>0.181731</td>
      <td>0.234064</td>
      <td>0.171499</td>
      <td>0.119430</td>
      <td>0.200905</td>
      <td>0.166378</td>
    </tr>
  </tbody>
</table>
</div>




```python
#Comparison analysis for ["Capomulin","Infubinol","Ketapril","Placebo" ]
#mean
Comp_metastatic_mean_Drug_effect=Metastatic_effect_mean_pivot_table[["Capomulin","Infubinol","Ketapril","Placebo" ]]
#sem
Comp_metastatic_sem_Drug_effect=Metastatic_effect_sem_pivot_table[["Capomulin","Infubinol","Ketapril","Placebo" ]]

```


```python
#Make a scatter plot with comparison analysis for ["Capomulin","Infubinol","Ketapril","Placebo" ] 
#based on the function above
colors = ['g','r','b','y']
markers=["^","s","d","o"]
color_index = 0
marker_index=0
for drug in Comp_mean_Drug_effect:
    scatterplot(x_data = Comp_metastatic_mean_Drug_effect.index, 
                y_data = Comp_metastatic_mean_Drug_effect[drug], error_data = Comp_metastatic_sem_Drug_effect [drug]
                ,color = colors[color_index],marker=markers[marker_index])
    color_index += 1
    marker_index +=1
    
xmin= Comp_metastatic_mean_Drug_effect.index.min()
xmax=Comp_metastatic_mean_Drug_effect.index.max()
ymin= pd.to_numeric(Comp_metastatic_mean_Drug_effect.min())
ymax=pd.to_numeric(Comp_metastatic_mean_Drug_effect.max())
#type(xmin)
plt.title("Metastatic Spread During Treatment")
plt.xlabel("Treatment Duration (Days)")
plt.ylabel("Met. Sites")
plt.legend(loc="best")
plt.grid(True, linestyle='dotted')
#plt.xticks(range(0,50, 5))
#plt.yticks(range(0.0,4.0,0.5)
plt.xlim(xmin, xmax)
plt.savefig("Fig2.png")
#plt.ylim(ymin, ymax)
plt.show()
```


![png](output_16_0.png)



```python
#Survival Rate
survival_rate= clinical_mouse_merge_df.groupby(["Drug","Timepoint"]).count()
survival_rate_df=pd.DataFrame(survival_rate['Mouse ID'])
survival_rate_df.head()

```




<div>
<style>
    .dataframe thead tr:only-child th {
        text-align: right;
    }

    .dataframe thead th {
        text-align: left;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th></th>
      <th>Mouse ID</th>
    </tr>
    <tr>
      <th>Drug</th>
      <th>Timepoint</th>
      <th></th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th rowspan="5" valign="top">Capomulin</th>
      <th>0</th>
      <td>25</td>
    </tr>
    <tr>
      <th>5</th>
      <td>25</td>
    </tr>
    <tr>
      <th>10</th>
      <td>25</td>
    </tr>
    <tr>
      <th>15</th>
      <td>24</td>
    </tr>
    <tr>
      <th>20</th>
      <td>23</td>
    </tr>
  </tbody>
</table>
</div>




```python
#Survival Rate pivot table
survival_rate_df_pivot_table=survival_rate_df.pivot_table('Mouse ID', ['Timepoint'], 'Drug')
survival_rate_df_pivot_table.head()
```




<div>
<style>
    .dataframe thead tr:only-child th {
        text-align: right;
    }

    .dataframe thead th {
        text-align: left;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>Drug</th>
      <th>Capomulin</th>
      <th>Ceftamin</th>
      <th>Infubinol</th>
      <th>Ketapril</th>
      <th>Naftisol</th>
      <th>Placebo</th>
      <th>Propriva</th>
      <th>Ramicane</th>
      <th>Stelasyn</th>
      <th>Zoniferol</th>
    </tr>
    <tr>
      <th>Timepoint</th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>25</td>
      <td>25</td>
      <td>25</td>
      <td>25</td>
      <td>25</td>
      <td>25</td>
      <td>26</td>
      <td>25</td>
      <td>26</td>
      <td>25</td>
    </tr>
    <tr>
      <th>5</th>
      <td>25</td>
      <td>21</td>
      <td>25</td>
      <td>23</td>
      <td>23</td>
      <td>24</td>
      <td>25</td>
      <td>25</td>
      <td>25</td>
      <td>24</td>
    </tr>
    <tr>
      <th>10</th>
      <td>25</td>
      <td>20</td>
      <td>21</td>
      <td>22</td>
      <td>21</td>
      <td>24</td>
      <td>23</td>
      <td>24</td>
      <td>23</td>
      <td>22</td>
    </tr>
    <tr>
      <th>15</th>
      <td>24</td>
      <td>19</td>
      <td>21</td>
      <td>19</td>
      <td>21</td>
      <td>20</td>
      <td>17</td>
      <td>24</td>
      <td>23</td>
      <td>21</td>
    </tr>
    <tr>
      <th>20</th>
      <td>23</td>
      <td>18</td>
      <td>20</td>
      <td>19</td>
      <td>20</td>
      <td>19</td>
      <td>17</td>
      <td>23</td>
      <td>21</td>
      <td>17</td>
    </tr>
  </tbody>
</table>
</div>




```python
#Number of mice at Timepoint "0" in each drug treatment
mouse_num_0=survival_rate_df_pivot_table.iloc[0,:]
#Pivot Table Percentage of survival rate per drug treatment
percentagem_survival_rate_df_pivot_table=(survival_rate_df_pivot_table/mouse_num_0)*100
#percentagem_survival_rate_df_pivot_table.head()
Comp_perc_survival_rate=percentagem_survival_rate_df_pivot_table[["Capomulin","Infubinol","Ketapril","Placebo" ]]
Comp_perc_survival_rate.head()
```




<div>
<style>
    .dataframe thead tr:only-child th {
        text-align: right;
    }

    .dataframe thead th {
        text-align: left;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>Drug</th>
      <th>Capomulin</th>
      <th>Infubinol</th>
      <th>Ketapril</th>
      <th>Placebo</th>
    </tr>
    <tr>
      <th>Timepoint</th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>100.0</td>
      <td>100.0</td>
      <td>100.0</td>
      <td>100.0</td>
    </tr>
    <tr>
      <th>5</th>
      <td>100.0</td>
      <td>100.0</td>
      <td>92.0</td>
      <td>96.0</td>
    </tr>
    <tr>
      <th>10</th>
      <td>100.0</td>
      <td>84.0</td>
      <td>88.0</td>
      <td>96.0</td>
    </tr>
    <tr>
      <th>15</th>
      <td>96.0</td>
      <td>84.0</td>
      <td>76.0</td>
      <td>80.0</td>
    </tr>
    <tr>
      <th>20</th>
      <td>92.0</td>
      <td>80.0</td>
      <td>76.0</td>
      <td>76.0</td>
    </tr>
  </tbody>
</table>
</div>




```python
#Make a scatter plot with comparison analysis for ["Capomulin","Infubinol","Ketapril","Placebo" ] 
#based on the function above
colors = ['g','r','b','y']
markers=["^","s","d","o"]
color_index = 0
marker_index=0
for drug in Comp_mean_Drug_effect:
    scatterplot(x_data =Comp_perc_survival_rate.index, 
                y_data = Comp_perc_survival_rate[drug], error_data=0
                ,color = colors[color_index],marker=markers[marker_index])
    color_index += 1
    marker_index +=1
    
xmin= Comp_perc_survival_rate.index.min()
xmax=Comp_perc_survival_rate.index.max()
#ymin= pd.to_numeric(Comp_metastatic_mean_Drug_effect.min())
#ymax=pd.to_numeric(Comp_metastatic_mean_Drug_effect.max())
#type(xmin)
plt.title("Survival During Treatment")
plt.xlabel("Time (Days)")
plt.ylabel("Survival Rate(%)")
plt.legend(loc="best")
plt.grid(True, linestyle='dotted')
#plt.xticks(range(0,50, 5))
#plt.yticks(range(0.0,4.0,0.5)
plt.xlim(xmin, xmax)
#plt.ylim(ymin, ymax)
plt.savefig("Fig3.png")
plt.show()
```


![png](output_20_0.png)



```python
#Creating a bar graph that compares the total % tumor volume change for each drug across the full 45 days.
Comp_mean_Drug_effect.set_index
#Tumor size at Timepoint = 0
tumor_size_timepoint_0=Comp_mean_Drug_effect.iloc[0,:]
#tumor_size_timepoint_45=Comp_mean_Drug_effect.loc[["45"],["Placebo"]]
tumor_size_timepoint_45=Comp_mean_Drug_effect.iloc[-1,:]
```


```python
variation_tumor_size_0_to_45=100*((tumor_size_timepoint_45/tumor_size_timepoint_0)-(tumor_size_timepoint_0/tumor_size_timepoint_0))
variation_tumor_size_0_to_45
```




    Drug
    Capomulin   -19.481481
    Infubinol    46.123457
    Ketapril     57.050505
    Placebo      51.272727
    dtype: float64




```python
x = variation_tumor_size_0_to_45.index
y = variation_tumor_size_0_to_45[x]
colors = []
for y in variation_tumor_size_0_to_45[x]:
    if y < 0:
        colors.append('g')
    else:
        colors.append('r')
colors
```




    ['g', 'r', 'r', 'r']




```python
plt.bar(x, variation_tumor_size_0_to_45, color=colors, align="center", width=0.5)
plt.grid(True, linestyle='dotted')
plt.title("Tumor Change Over 45 Days Treatment")
plt.ylabel("%Tumor Volume Change")


plt.savefig("Fig4.png")
plt.show()
```


![png](output_24_0.png)



```python
#Have no idea how to add the percentages to the plot, have to ask!


```
