# FM1-Project
Fluid Mechanics 1 Project

ابتدا فایل Catalog_data.rar را اکسترکت کنید. هگام اجرای کد pump_detector.m (تشخیص پمپ مناسب)، اطمینان حاصل کنید که پوشه مربوط به داده‌های کاتالوگ (مانند فایل‌های عملکردی پمپ‌ها، قطر پروانه، توان و بازده) با عنوان Catalog_data در مسیر فعلی (Current Directory) قرار دارد یا در مسیر addpath شده باشد. در غیر این صورت، کد قادر به خواندن فایل‌های .xls یا .xlsx نخواهد بود.
First extract the Catalog.rar file. When running the pump_detector.m script (Pump Selector), make sure that the folder Catalog_data (containing pump performance curves, diameter, power, and efficiency files) is in the current MATLAB directory or properly added to the path. Otherwise, the code won't be able to access the required .xls/.xlsx files.

ه
Data File Naming Convention       فرمت نامگذاری فایل های دیتا اکسل .
| Format    | Description                          |
| --------- | ------------------------------------ |
| `%-%`     | Flow–Head performance (Q–H)          |
| `d-%-%`   | List of impeller diameters available |
| `%-%-%`   | Head curve at specific diameter      |
| `%-%---%` | Power data at specific diameter      |
| `e-%-%`   | Available efficiency levels          |
| `%-%--%`  | Efficiency curve at a fixed value    |

