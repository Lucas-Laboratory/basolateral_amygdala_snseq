# Drop-in patch instructions

**Tested in `scvi-tools v1.3.3`.**

*This patch exposes several `pytorch` loss functions for logging and plotting.*

## Instructions

To patch scvi_tools, replace the default `_base_model.py` file with the file contained in this directory.

```path
scvi_tools_hotpatch/
   └──src/
      └──scvi/
         └──base/
            └──model/
               └──base_model.py
```
*File path follows `scvi-tools` architecture*


**This repository includes code derived from: Yosef Lab, Weizmann Institute of Science**
**BSD 3-Clause License (see `pipelines/scvi_tools_cli_studio/module/scvi_tools_drop_in/LICENSE.bsd3-clause`)**