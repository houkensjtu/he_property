## he_property
> Helium-4 property calculation experiments.

NIST发布的Helium物性计算程序，以及Regen中都用到了Helium物性计算的Fortran程序。这些程序的依据来自于NIST的研究者Robert D. McCarty和Vincent D. Arp
发布的研究报告以及论文"A New Wide Range Equation of State for Helium"。最近，在2003年由一个Texas A&M大学的PhD学生重新更新了这些物性计算的方法，
在论文中他注明MaCarty和Arp都参与了指导他的论文工作，所以文章的可信程度较高。   
由于NIST最初的研究报告暂时找不到入手途径，以下是我学习这篇论文的一些笔记：    
### 1. Introduction
#### 1.1 Helium properties
- Helium发现于1868年，Helium这个词语来自于希腊语，意为太阳。
- Helium在整个宇宙中的比例仅次于氢气，但是在地球大气中极为稀有。Helium有两种同位体，He4和He3，这篇论文仅关注于He4的物性计算。
- Helium具有**所有气体中最低的沸点(4.2K在一个大气压）**，并在接近0K时都不会固化。氢气的沸点为20.7K，氖为27.1K。**这种性质使Helium成为了理想的制冷媒体。**
- Helium的第一个工业应用是充填气球。氢气可以提供更大的升力，但是由于Helium更加稳定不会爆炸，使之成为理想的气体选择。
- Helium在所有气体中**具有较高的热传导率，同时极易扩散**，使之在光纤的制造中被应用为环境气体。
- Helium的低沸点使其在各种低温应用中有特殊地位。低温应用包括制冷机，超导磁体，加速器，MRI以及SQUID等医疗仪器。**低温应用是目前消费Helium最多的工业应用**。
- Helium也被用来purge和pressurize氢气推进器（火箭及防御用途），因为只有Helium的沸点低于氢气，使用其他气体会导致冻结并损坏引擎。
- Helium同时具有特殊的相变化性质，在2.2K以下会转移成superfluid。因此在学术上也具有特殊的研究价值。

#### 1.2 Supply and usage
- Helium在一战以前全部由政府出资挖掘。在一战中发现Helium的重要性以及不足之后，1960年的法案通过允许私人企业开采挖掘Helium。
- 在美国国内，开采Helium的主要field在Texas州附近；其形式是通过从天然气中提炼。全世界在美国之外，只有2处可以开采Helium的地方，所以全世界的Helium应用
高度依赖于美国。
- 调查显示美国国内的Helium储量开始下降，同时在Russia，非洲以及澳大利亚有新的小型开采源发现。全球的helium产量开始不那么依赖于美国。
- 美国拥有世界上最庞大先进的超导磁体研究/应用的重要原因，是美国拥有廉价的国内Helium开采途径。半导体制造过程中也需要helium。

#### 1.3 IMPORTANCE OF ACCURATE THERMODYNAMIC PROPERTIES FOR HELIUM
- Helium的特殊物性使其在很多应用中**具有无法替代的地位**；另一方面，Helium被使用后无法再生，在全球面临耗竭的危险。美国在2012年颁布行动令“Helium Stewardship Act of 2012”，鼓励研究提炼，再生Helium的新方法。
- 精确的物性计算是高效利用Helium，设计Helium机械所必须的情报。
- Helium物性最初由实验室的测量数据产生，这些数据只能用来查找，却不便于利用于计算机程序。之后，很多研究者发布了拟合这些数据的拟合方法，这些方法
大多局限于某些特定物性，或是特定的温度压力领域，不够具有广泛性。
- 

