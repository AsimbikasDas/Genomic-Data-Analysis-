<h1><b></b>GenomeLREP</b></h1>

<h2><b>Description:</b></h2>

GenomeLRep is a web-based platform designed to process and analyze genomic datasets through advanced regularization techniques. The project involved developing an interactive website that allows users to upload and work with complex genomic data efficiently.

We implemented various regularization algorithms such as Lasso, Ridge Regression, Elastic Net, and Principal Component Analysis (PCA) to clean, compress, and prepare the data for deeper statistical and machine learning analysis. These techniques helped to reduce overfitting, manage multicollinearity, and extract meaningful patterns from high-dimensional genetic data.

Once the data was regularized, it was further utilized for regression modeling to identify relationships between genetic variables and target outcomes. Additionally, the platform includes features for graphical visualization, enabling users to explore data trends and model results through dynamic plots and charts.

The goal of GenomeLRep is to simplify the process of genomic data preparation and modeling, making powerful bioinformatics tools more accessible to researchers, data scientists, and healthcare professionals.

<h2>🚀 Getting Started</h2>
<b>Most important thing to note is that this project was bulit by our team from scratch and also tested in Pycharm(Educational License) but you can use any other IDEs like Visual Studio code or other softwares if any, whatever you are most comfortable with.</b>

Follow these steps to set up and run the project in an isolated Python virtual environment:

<h3>1. Clone the Repository</h3>
Open your terminal or command prompt and run:

```
git clone https://github.com/Umamaheshwarmydam/GenomeLREP1.git
cd GenomeLREP1
```

<h3>2. Set Up a Virtual Environment</h3>
Creating a virtual environment ensures that project dependencies are isolated from your global Python installation.​
<li>On Windows:</li>

```
python -m venv venv
venv\Scripts\activate
```

<li>On macOS/Linux:</li>

```
python3 -m venv venv
source venv/bin/activate
```

Once activated, your terminal prompt will change to indicate that you're working within the virtual environment.​

<h3>3. Install Dependencies</h3>
With the virtual environment activated, install the required packages:​

```
pip install -r requirements.txt
```

This command installs all the dependencies specified in the requirements.txt file.​

<h3>4. Run the Flask Application</h3>
Start the Flask development server:

```
python app.py
```

By default, the application will be accessible at http://127.0.0.1:5000/ in your web browser.

```
GenomeLREP1/
├── app.py
├── requirements.txt
├── static/
├── templates/
├── README.md
└── .idea/
```

<h2><b>A Guide to use our website:</b></h2>

![WhatsApp Image 2025-04-24 at 14 17 45_5097e31e](https://github.com/user-attachments/assets/f07bbbeb-5301-4b0d-a872-54affd1bed2a)

<pre>
<h3><b>The above image is the home/landing page of our website.</b></h3>

</pre>

[![WhatsApp Image 2025-04-24 at 14 19 07_a86a3b59](https://github.com/user-attachments/assets/090c9183-3d2c-406f-9803-f76422f3aa4b)](https://ncase.me/loopy/v1.1/?embed=1&data=[[[1,1228,-34,1,%22User%2520input%2520%22,4],[2,1188,229,1,%22Website%2520%22,5],[3,856,365,0.5,%22Backend%2520%22,2],[5,427,607,0.5,%22ML%2520library%22,4],[12,387,931,0.5,%22LASSO%22,4],[13,784,896,0.5,%22Ridge%22,4],[14,1138,793,0.5,%22Elastic%2520net%22,4],[15,1130,550,0.5,%22JSON%2520FORMAT%22,3],[16,519,25,0.5,%22Front%2520end%22,1],[21,63,715,0.5,%22PCA%22,0]],[[2,1,44,-1,0],[1,2,40,1,0],[2,3,-43,1,0],[3,2,-92,-1,0],[5,3,141,1,0],[12,5,70,1,0],[13,5,-74,1,0],[14,5,-87,1,0],[5,15,68,1,0],[15,2,-115,1,0],[2,16,32,1,0],[16,2,53,-1,0],[3,5,-33,-1,0],[5,12,20,-1,0],[5,13,-15,-1,0],[5,14,44,-1,0],[5,21,30,-1,0],[21,5,94,1,0]],[],21%5D)

<pre>
<h3><b>The above image is the flow diagram of how our interface system works, give a tap on it then you will be 
redirected to another website there tap on frontend node and that's all you can visualise how our 
interface works.</b></h3>

</pre>

![WhatsApp Image 2025-04-24 at 14 21 03_f98e70f4](https://github.com/user-attachments/assets/18036b98-159f-442b-9454-9778bc668841)

<pre>
<h3><b>The above image represents our services page where you can see file uploading feature and selecting 
methods which you want to be applied on your dataset.</b></h3>

</pre>

![WhatsApp Image 2025-04-24 at 14 23 03_8bb41e7f](https://github.com/user-attachments/assets/cc5386fb-7e85-4b1a-803a-7780d8e65459)

<pre>
<h3><b>The above shows a file uploaded and when selected all methods we will get to see hyperparameters of all 
methods if you are aware of what parameters you have to select you can do so. Otherwise you can keep it 
default and press blue colored Analyze data button.</b></h3>

</pre>


![WhatsApp Image 2025-04-24 at 14 24 25_d11f4716](https://github.com/user-attachments/assets/5989ad0a-f80b-4463-ad79-f9405a58186c)

<pre>
<h3><b>In the above we can see a separate button to download Regression methods and PCA values table as well as 
another button for downloading pca plot image for further research this downloadable files can be used in 
bioinformatics. After this one can run new analysis by clicking right most button and do the same process 
as described early</b></h3>
  
</pre>

![WhatsApp Image 2025-04-24 at 14 25 49_2da3ca5b](https://github.com/user-attachments/assets/6e0ee95e-028a-42f6-8dc0-bd31717f116b)

<pre>
<h3><b>The above image gives information about our team members who contributed to this project.</b></h3>
  
</pre>

<h1><b>📌 Notes</b></h1>
<li>Ensure that you have Python installed on your system before proceeding.</li>

<li>Using a virtual environment is a best practice to manage project-specific dependencies and avoid conflicts.</li>
