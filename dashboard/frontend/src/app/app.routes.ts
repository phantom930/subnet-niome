import { Routes } from '@angular/router';

export const routes: Routes = [
  { path: '', pathMatch: 'full', redirectTo: 'tasks' },
  {
    path: 'tasks',
    title: 'Tasks',
    loadComponent: () => import('./pages/tasks/tasks').then((m) => m.Tasks),
  },
  {
    path: 'mining',
    title: 'Mining process',
    loadComponent: () => import('./pages/mining/mining').then((m) => m.Mining),
  },
  { path: '**', redirectTo: 'tasks' },
];
